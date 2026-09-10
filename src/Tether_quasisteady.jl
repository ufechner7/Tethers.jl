# A quasi-steady tether model: solves for tether shape and forces given the ground-station
# orientation/tension and the kite's position and velocity.
#
# This is a submodule of `Tethers`, so that its names do not end up in `Main`: several
# examples (e.g. `Tether_08.jl`) define their own top-level `Settings` struct, and
# `runtests.jl` includes all of them into the same `Main`, which a top-level `Settings`
# here would collide with. Use it with `using Tethers.QuasiSteady`.
module QuasiSteady

using LinearAlgebra, StaticArrays, ADTypes, NonlinearSolve, MAT, Parameters#, QuadGK

export StaticSettings, Tether, init!, step!, clear!, elevation, azimuth, tension,
       n_nodes, get_initial_conditions, get_analytic_catenary

include(joinpath(@__DIR__, "qsm_conventions.jl"))

const MVec3 = MVector{3, Float64}
const SVec3 = SVector{3, Float64}

# Solvers used by [`simulate_tether`](@ref). This problem needs globalization: the initial
# guess for the tension is easily orders of magnitude off, and an undamped Newton step then
# jumps clean out of the physically meaningful region - straight to `Tn` = 0, where the
# residual flattens out at the shape of a tether hanging limp from the ground station and
# there is no gradient left to come back on. The radius is capped because the unknowns are
# O(1) by construction (see `scaled_res`), so a step longer than 2 is never useful, while
# the default cap of `max(norm(f(u0)), maximum(u0) - minimum(u0))` lets a single step cross
# ten decades of tension. `AutoForwardDiff` gives an exact 3x3 Jacobian in one evaluation,
# where `AutoFiniteDiff` needs four inexact ones.
const DEFAULT_SOLVER = TrustRegion(autodiff = AutoForwardDiff(),
                                   max_trust_radius = 2.0, initial_trust_radius = 1.0)
# Only used if `DEFAULT_SOLVER` fails. The linear parameterization of the tension is worse
# conditioned, but it fails in different places, which is the point of a fallback.
const FALLBACK_SOLVER = TrustRegion(autodiff = AutoForwardDiff())

@with_kw mutable struct StaticSettings @deftype Float64
    "number of tether segments; the model stores `segments - 1` nodes, see [`n_nodes`](@ref)"
    segments::Int64 = 8
    "initial elevation angle β                                     [deg]"
    elevation = 70.0
    "initial wind-frame azimuth angle φ                            [deg]"
    azimuth = 0.0
    "initial unstretched tether length                               [m]"
    l_tether = 50.0
    "initial tether slack, l_tether = (1 + slack) * kite distance"
    slack = 0.05
    "density of air                                              [kg/m³]"
    rho = 1.225
    "gravitational acceleration, W frame                          [m/s²]"
    g_earth::MVector{3, Float64} = [0.0, 0.0, -9.81]
    "drag coefficient of the tether"
    cd_tether = 0.958
    "diameter of the tether                                          [mm]"
    d_tether = 4
    "density of the tether (Dyneema)                              [kg/m³]"
    rho_tether = 724
    "unit spring constant of the tether (= E*A)                       [N]"
    c_spring = 614600
    "the nonlinear solver used by init!/step!"
    alg::Any = DEFAULT_SOLVER
end

"""
    StaticSettings

Physical, structural and solver parameters of the quasi-steady tether model, plus the
initial condition used by [`init!`](@ref). Nothing in it changes while a simulation
runs - see [`Tether`](@ref) for the state that does.

# Fields
  - segments::Int64: number of tether segments; `segments - 1` nodes are stored, see [`n_nodes`](@ref)
  - elevation::Float64: initial elevation angle β [deg]
  - azimuth::Float64: initial wind-frame azimuth angle φ [deg]
  - l_tether::Float64: initial unstretched tether length [m]
  - slack::Float64: initial tether slack, `l_tether = (1 + slack) * kite distance`
  - rho::Float64: density of air [kg/m³]
  - g_earth::MVector{Float64}: gravitational acceleration [m/s²]
  - cd_tether::Float64: drag coefficient of the tether
  - d_tether::Float64: diameter of the tether [mm]
  - rho_tether::Float64: density of the tether (Dyneema) [kg/m³]
  - c_spring::Float64: unit spring constant [N] (= `E*A`)
  - alg: the nonlinear solver used by [`init!`](@ref)/[`step!`](@ref), defaults to `DEFAULT_SOLVER`
"""
StaticSettings

"""
    n_nodes(segments)

Number of tether nodes the model stores for a tether of `segments` segments.

`segments` counts segments, as in the dynamic mass-spring model, so a tether of `segments`
segments has `segments + 1` points. Neither of the two outer ones is stored: the ground
station sits at the origin and the kite attachment point is returned separately as `p0`,
which leaves `segments - 1` nodes in between.
"""
@inline n_nodes(segments::Integer) = segments - 1

"""
    Tether

Mutable state of one quasi-steady tether simulation: the settings, the persistent
solver state, the boundary conditions of the last [`step!`](@ref) and its results. All
vectors are in the W (wind) reference frame. Buffers are sized from `set.segments` at
construction and reused by every `step!`, so a stepping loop stays allocation free.

Construct with `Tether(se::StaticSettings)`, then call [`init!`](@ref) before the first
[`step!`](@ref).

# Fields
  - set::StaticSettings: settings, see [`StaticSettings`](@ref)
  - state_vec::MVector{3, Float64}: (β [rad], φ [rad], Tn [N]) at the ground station
  - kite_pos::MVector{3, Float64}: kite position of the last step [m]
  - kite_vel::MVector{3, Float64}: kite velocity of the last step [m/s]
  - wind_vel::Matrix{Float64}: wind velocity per segment of the last step [m/s]
  - tether_length::Float64: unstretched tether length of the last step [m]
  - tether_pos::Matrix{Float64}: node coordinates of the last step [m]
  - force_gnd::Float64: tension at the ground station [N]
  - force_kite::MVector{3, Float64}: force on the tether at the kite attachment point [N]
  - p0::MVector{3, Float64}: kite-tether attachment point [m]
"""
@with_kw mutable struct Tether
    set::StaticSettings = StaticSettings()
    # persistent state, updated by init! and step!
    state_vec::MVector{3, Float64} = zeros(MVector{3})   # (β [rad], φ [rad], Tn [N])
    # boundary conditions of the last step
    kite_pos::MVector{3, Float64} = zeros(MVector{3})
    kite_vel::MVector{3, Float64} = zeros(MVector{3})
    wind_vel::Matrix{Float64} = zeros(3, n_nodes(set.segments))
    tether_length::Float64 = 0.0
    # results of the last step
    tether_pos::Matrix{Float64} = zeros(3, n_nodes(set.segments))  # node coordinates
    force_gnd::Float64 = 0.0                              # tension at ground station [N]
    force_kite::MVector{3, Float64} = zeros(MVector{3})   # force on the tether end
    p0::MVector{3, Float64} = zeros(MVector{3})           # kite-tether attachment point
end

"""
    Tether(se::StaticSettings)

Construct a `Tether` whose buffers are sized for `se.segments`. The state is left at
zero until [`init!`](@ref) runs.
"""
Tether(se::StaticSettings) = Tether(set=se)

"""
    elevation(te::Tether)

Elevation angle β [rad] at the ground station, from `te.state_vec`.
"""
elevation(te::Tether) = te.state_vec[1]

"""
    azimuth(te::Tether)

Wind-frame azimuth angle φ [rad] at the ground station, from `te.state_vec`.
"""
azimuth(te::Tether) = te.state_vec[2]

"""
    tension(te::Tether)

Tension Tn [N] at the ground station, from `te.state_vec`.
"""
tension(te::Tether) = te.state_vec[3]

"""
    clear!(te::Tether)

Reset `te`'s persistent state and result buffers to zero, resizing them if
`te.set.segments` has changed since construction. Does not touch `te.set`. Called by
[`init!`](@ref); also useful on its own to restart a simulation with the same `Tether`.

# Returns
- te::Tether, cleared
"""
function clear!(te::Tether)
    n = n_nodes(te.set.segments)
    te.state_vec .= 0.0
    te.kite_pos .= 0.0
    te.kite_vel .= 0.0
    size(te.wind_vel, 2) == n ? (te.wind_vel .= 0.0) : (te.wind_vel = zeros(3, n))
    te.tether_length = 0.0
    size(te.tether_pos, 2) == n ? (te.tether_pos .= 0.0) : (te.tether_pos = zeros(3, n))
    te.force_gnd = 0.0
    te.force_kite .= 0.0
    te.p0 .= 0.0
    te
end

"""
    check_wind_vel(wind_vel, segments)

Validate that `wind_vel` is a `(3, n_nodes(segments))` matrix - one column per stored node,
see [`n_nodes`](@ref) - as required by [`init!`](@ref) and [`step!`](@ref). `segments`
always comes from `te.set.segments` - it is never inferred from `size(wind_vel, 2)`,
unlike the internal [`simulate_tether`](@ref).
"""
function check_wind_vel(wind_vel, segments)
    size(wind_vel, 1) == 3 ||
        throw(ArgumentError("wind_vel must have 3 rows, got $(size(wind_vel, 1))"))
    size(wind_vel, 2) == n_nodes(segments) ||
        throw(ArgumentError("wind_vel must have n_nodes(te.set.segments) = " *
                             "$(n_nodes(segments)) columns, got $(size(wind_vel, 2))"))
    nothing
end

"""
    init!(te::Tether; prn = false)

Initialize `te`: derive the initial kite position from `te.set.elevation`,
`te.set.azimuth` and `te.set.l_tether`, solve the catenary equation for an initial
guess of `te.state_vec`, then run one [`step!`](@ref) so that `te` is left in a
consistent, solved state. Takes no state arguments - everything comes from `te.set`.

# Keyword arguments
- prn: print the solver statistics of the final `step!`

# Returns
- te::Tether, initialized and solved
"""
function init!(te::Tether; prn = false)
    clear!(te)
    se = te.set
    check_wind_vel(te.wind_vel, se.segments)
    β0, φ0 = deg2rad(se.elevation), deg2rad(se.azimuth)
    kite_distance = se.l_tether / (1 + se.slack)
    kite_pos = MVector{3}(kite_distance * cos(β0) * cos(φ0),
                           kite_distance * cos(β0) * sin(φ0),
                           kite_distance * sin(β0))

    state_vec, _, _, _, _, _ = init_quasisteady(kite_pos, se.l_tether; kite_vel=te.kite_vel,
                                                 segments=se.segments, wind_vel=te.wind_vel,
                                                 settings=se)
    te.state_vec .= state_vec

    step!(te, kite_pos, te.kite_vel; tether_length=se.l_tether, wind_vel=te.wind_vel, prn)
end

"""
    step!(te::Tether, kite_pos, kite_vel; tether_length = nothing, wind_vel = nothing, prn = false)

Move the loose end of the tether to `kite_pos`/`kite_vel` and re-solve for the tether
shape and forces, using `te.state_vec` as the initial guess. Writes `state_vec`,
`kite_pos`, `kite_vel`, `wind_vel`, `tether_length`, `tether_pos`, `force_gnd`,
`force_kite` and `p0` into `te`.

# Arguments
- te::Tether: the tether, see [`Tether`](@ref)
- kite_pos::MVector{3, Float64}: kite position in the W frame [m]
- kite_vel::MVector{3, Float64}: kite velocity in the W frame [m/s]

# Keyword arguments
- tether_length: unstretched tether length [m]; defaults to
  `(1 + te.set.slack) * norm(kite_pos)`
- wind_vel: `(3, n_nodes(te.set.segments))` matrix, wind velocity per node [m/s]; defaults to
  `te.wind_vel`
- prn: print the solver statistics

# Returns
- te::Tether, with all fields above updated
"""
function step!(te::Tether, kite_pos, kite_vel; tether_length=nothing, wind_vel=nothing,
               prn=false)
    se = te.set
    _tether_length = tether_length === nothing ? (1 + se.slack) * norm(kite_pos) : tether_length
    _wind_vel = wind_vel === nothing ? te.wind_vel : wind_vel
    check_wind_vel(_wind_vel, se.segments)

    state_vec, _, force_gnd, force_kite, p0 = simulate_tether(
        te.state_vec, kite_pos, kite_vel, _wind_vel, _tether_length, se;
        prn, alg=se.alg, tether_pos=te.tether_pos)

    te.state_vec .= state_vec
    te.kite_pos .= kite_pos
    te.kite_vel .= kite_vel
    te.wind_vel = _wind_vel
    te.tether_length = _tether_length
    # te.tether_pos was written in place by simulate_tether via the `tether_pos` keyword
    te.force_gnd = force_gnd
    te.force_kite .= force_kite
    te.p0 .= p0
    te
end

"""
    simulate_tether(state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings)

Function to determine the tether shape and forces, based on a quasi-steady model.

# Arguments
- state_vec::MVector{3, Float64}: state vector (beta [rad], phi [rad], Tn [N]);
  tether orientation and tension at ground station
- kite_pos::MVector{3, Float64}: kite position vector in wind reference frame
- kite_vel::MVector{3, Float64}: kite velocity vector in wind reference frame
- wind_vel:: (3, n_nodes(segments)) MMatrix{Float64} wind velocity vector in wind reference frame, one column per node
- tether_length: tether length
- settings:: StaticSettings struct containing environmental and tether parameters: see [`StaticSettings`](@ref)

# Keyword arguments
- prn: print the solver statistics
- alg: the nonlinear solver, defaults to `DEFAULT_SOLVER`
- tether_pos: `(3, n_nodes(segments))` matrix that receives the node positions, or `nothing` to
  allocate a fresh one (the default). Passing a pre-allocated buffer, as [`step!`](@ref)
  does, avoids that allocation in a stepping loop.

# Returns
- state_vec::MVector{3, Float64}: state vector (beta [rad], phi [rad], Tn [N]);
  tether orientation and tension at ground station
- tether_pos::Matrix{Float64}: x,y,z - coordinates of the tether nodes
- force_gnd::Float64: Line tension at the ground station
- force_kite::MVector{3, Float64}: force from the kite to the end of tether
- p0::MVector{3, Float64}:  x,y,z - coordinates of the kite-tether attachment
"""
function simulate_tether(state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings;
                         prn=false, alg=DEFAULT_SOLVER, tether_pos=nothing)
    segments = size(wind_vel, 2) + 1   # wind_vel has one column per node, see `n_nodes`
    # The tension is solved for as log(Tn/tension_scale), see `scaled_res`. `tension_scale`
    # is the initial guess itself, so the third unknown simply starts at zero.
    tension_scale = state_vec[3] > 0 ? Float64(state_vec[3]) : 2e-4 * settings.c_spring
    param = (kite_pos=SVec3(kite_pos), kite_vel=SVec3(kite_vel), wind_vel=wind_vel,
             tether_length=tether_length, settings=settings, segments=segments,
             tension_scale=tension_scale)

    # Out-of-place residual on a static vector, which keeps the solver allocation free.
    u0 = SVec3(state_vec[1], state_vec[2], 0.0)
    prob = NonlinearProblem{false}(scaled_res, u0, param)
    sol = solve(prob, alg)
    tol = 1e-6 * tether_length
    if converged(sol, tol)
        tension = tension_scale * exp(sol.u[3])
    else
        sol = solve(NonlinearProblem{false}(lin_res, SVec3(u0[1], u0[2], 1.0), param),
                    FALLBACK_SOLVER)
        tension = tension_scale * sol.u[3]
        converged(sol, tol) ||
            @warn "simulate_tether did not converge" sol.retcode norm(sol.resid)
    end
    if prn
        nsteps = sol.stats === nothing ? "n/a" : sol.stats.nsteps
        println("Iterations: ", nsteps, ", retcode: ", sol.retcode, ", |res|: ", norm(sol.resid))
    end

    state_vec = MVector(sol.u[1], sol.u[2], tension)
    # Re-run the model at the solution, this time storing the node positions.
    tether_pos_buf = tether_pos === nothing ? Matrix{Float64}(undef, 3, n_nodes(segments)) : tether_pos
    _, force_kite, p0 = tether_shape(state_vec[1], state_vec[2], state_vec[3], param, tether_pos_buf)

    force_gnd = state_vec[3]
    state_vec, tether_pos_buf, force_gnd, MVector(force_kite), MVector(p0)
end

"""
    scaled_res(u, param)

Residual of [`tether_shape`](@ref) as seen by the nonlinear solver: `u` is
`(beta, phi, log(Tn/tension_scale))`, see [`simulate_tether`](@ref).

The tension is solved for on a logarithmic scale because it spans decades. A taut tether
pulls with 1e5 N, a tether long enough to sag between kite and ground station with a few N,
so an initial guess can easily be a factor 1e5 off, and a solver that walks that distance
linearly needs an iteration per decade. On top of that the residual is far less sensitive
to the tension than to the two angles - a stiff tether hardly stretches - which leaves the
Jacobian nearly singular in the tension direction; `d(res)/d(log Tn)` is `Tn * d(res)/dTn`,
which is of the same order as the two angle columns. Positivity of the tension comes free.
"""
scaled_res(u, param) = first(tether_shape(u[1], u[2], param.tension_scale * exp(u[3]), param, nothing))

"""
    lin_res(u, param)

As [`scaled_res`](@ref), but with `u[3]` the tension itself in units of `tension_scale`.
Used only by the fallback solve in [`simulate_tether`](@ref).
"""
lin_res(u, param) = first(tether_shape(u[1], u[2], u[3] * param.tension_scale, param, nothing))

"""
    converged(sol, tol)

Whether a nonlinear solution actually solved the problem. Written so that a `NaN` residual
counts as "not converged" rather than slipping through a `>` comparison.
"""
converged(sol, tol) = sol.retcode == ReturnCode.Success && norm(sol.resid) < tol

"""
    node_kinematics(ω, p_unit, v_parallel, pos)

Velocity and acceleration of a tether node at `pos`, assuming the tether rotates rigidly
with the kite: `ω` is the angular velocity of the kite position vector, and `v_parallel` the
radial component of the kite velocity along the unit vector `p_unit`.
"""
@inline function node_kinematics(ω, p_unit, v_parallel, pos)
    vel = v_parallel * p_unit + cross(ω, pos)
    acc = cross(ω, cross(ω, pos))
    return vel, acc
end

"""
    segment_drag(v_app, dir, drag_coeff)

Drag force on one tether segment with unit direction `dir` and apparent wind velocity
`v_app`: `drag_coeff * |v_n| * v_n`, with `v_n` the component of `v_app` normal to the
segment. Below 1 mm/s of apparent wind the drag is zero by definition, which also keeps the
normal direction well defined.
"""
@inline function segment_drag(v_app, dir, drag_coeff)
    if abs(v_app[1]) < 1e-3 && abs(v_app[2]) < 1e-3 && abs(v_app[3]) < 1e-3
        return zero(v_app)
    end
    v_n = v_app - dot(v_app, dir) * dir
    n2 = dot(v_n, v_n)
    n2 < 1e-24 && return zero(v_app)   # apparent wind aligned with the segment
    return (drag_coeff * sqrt(n2)) * v_n
end

@inline function set_col!(m, i, v)
    @inbounds m[1, i], m[2, i], m[3, i] = v[1], v[2], v[3]
    return nothing
end

@inline wind_col(w, i) = @inbounds SVector(w[1, i], w[2, i], w[3, i])

"""
    tether_shape(β, φ, Tn, param, pj)

Integrate the quasi-steady tether from the ground station up to the kite and return the
gap between the kite and the end of the tether.

The tether is walked one segment at a time, so only the force, drag, velocity and
acceleration of the *current* node are needed; keeping them in `SVector`s instead of in
`(3, n_nodes(segments))` buffers makes the whole integration allocation free and lets ForwardDiff
run straight through it.

# Arguments
- β, φ, Tn: elevation [rad], wind-frame azimuth [rad] and tension [N] at the ground station
- param: named tuple with `kite_pos`, `kite_vel`, `wind_vel`, `tether_length`, `settings`
  and `segments`, see [`simulate_tether`](@ref)
- pj: `(3, n_nodes(segments))` matrix that receives the node positions, or `nothing` to skip
  them. Node `n_nodes(segments)` is the one closest to the ground station, node 1 the last one before the
  kite attachment point `p0`.

# Returns
- res: difference between the kite position and the end of the tether `p0`
- T0: force from the kite on the end of the tether
- p0: x,y,z - coordinates of the kite-tether attachment
"""
function tether_shape(β, φ, Tn, param, pj)
    (; kite_pos, kite_vel, wind_vel, tether_length, settings, segments) = param
    g = abs(settings.g_earth[3])
    nn = n_nodes(segments)          # stored nodes; there is one segment more than that
    Ls = tether_length / segments
    # the reference area of a segment is its length times its diameter, in m²; `d_tether` is
    # in mm, so it needs the same /1000 as the cross section below
    drag_coeff = -0.5 * settings.rho * Ls * (settings.d_tether/1000) * settings.cd_tether
    A = π/4 * (settings.d_tether/1000)^2
    mj = settings.rho_tether * Ls * A
    EA = settings.c_spring          # the model's E is c_spring/A, so E*A is c_spring again

    # Precompute common values
    sinβ, cosβ = sin(β), cos(β)
    sinφ, cosφ = sin(φ), cos(φ)
    kite_p = SVec3(kite_pos)
    kite_v = SVec3(kite_vel)
    norm_p = norm(kite_p)
    p_unit = kite_p / norm_p
    v_parallel = dot(kite_v, p_unit)
    ω = cross(kite_p / norm_p^2, kite_v)

    # First element: the segment leaving the ground station, ending in node `nn`
    dir = SVector(cosβ*cosφ, cosβ*sinφ, sinβ)   # cos(elevation)cos(azimuth), ...
    FT = SVector(Tn*cosβ*cosφ, Tn*cosβ*sinφ, Tn*sinβ)
    pos = Ls * dir
    pj === nothing || set_col!(pj, nn, pos)
    vel, acc = node_kinematics(ω, p_unit, v_parallel, pos)
    Fd = segment_drag(vel - wind_col(wind_vel, nn), dir, drag_coeff)

    # Process the other segments, walking up towards the kite
    @inbounds for ii in nn:-1:2
        # Tension force: the node below carries 1.5 segment masses (the model lumps the
        # half segment at the ground station onto it), all others carry one.
        mj_total = ii == nn ? 1.5mj : mj
        FT = mj_total * acc + FT - Fd + SVector(0.0, 0.0, mj_total * g)

        # Position of the next node, the segment being stretched by its own tension
        ft_norm = norm(FT)
        l_i_1 = (ft_norm/EA + 1) * Ls
        pos_next = pos + l_i_1 * (FT / ft_norm)

        # The drag of this segment uses the velocity of the node below it, so it has to be
        # taken before `vel` is advanced.
        v_app = vel - wind_col(wind_vel, ii)
        vel, acc = node_kinematics(ω, p_unit, v_parallel, pos_next)
        seg = pos_next - pos
        Fd = segment_drag(v_app, seg / norm(seg), drag_coeff)

        pos = pos_next
        pj === nothing || set_col!(pj, ii-1, pos)
    end

    # Final ground connection calculations
    T0 = 1.5mj * acc + FT - Fd + SVector(0.0, 0.0, 1.5mj * g)
    T0_norm = norm(T0)
    l_i_1 = (T0_norm/EA + 1) * Ls
    p0 = pos + l_i_1 * (T0 / T0_norm)

    return kite_p - p0, T0, p0
end

"""
    res!(res, state_vec, param)

Calculates difference between tether end and kite given tether ground segment orientation 
and magnitude. Thin, mutating wrapper around [`tether_shape`](@ref), kept for callers that
work with the in-place `(res, state_vec, param)` signature.

# Arguments
- res::Vector{Float64} difference between tether end and kite segment
- state_vec::MVector{3, Float64} state vector (beta [rad], phi [rad], Tn [N]);
  tether orientation and tension at ground station
- par:: 8-elements tuple:
    - kite_pos::MVector{3, Float64} kite position vector in wind reference frame
    - kite_vel::MVector{3, Float64} kite velocity vector in wind reference frame
    - wind_vel::MMatrix{Float64} wind velocity vector in wind reference frame for each segment of the tether
    - tether_length: tether length
    - settings:: StaticSettings struct containing environmental and tether parameters: see [`StaticSettings`](@ref)
    - buffers:: (5, ) Vector{Matrix{Float64}}  Vector of (3, n_nodes(segments)) Matrix{Float64} empty matrices;
      only `buffers[3]` is used, it receives the node positions
    - segments:: number of tether segments; `segments - 1` nodes are stored
    - return_result:: Boolean to determine use for in-place optimization or for calculating returns

# Returns (if return_result==true)
- res::Vector{Float64} difference between tether end and kite segment
- T0::MVector{3, Float64} force from the kite to the end of tether
- pj:: (3, n_nodes(segments)) Matrix{Float64} x,y,z - coordinates of the tether nodes
- p0::MVector{3, Float64}  x,y,z - coordinates of the kite-tether attachment

# Example usage
```julia
state_vec = rand(3,)
kite_pos = [100, 100, 300] 
kite_vel = [0, 0, 0]
wind_vel = rand(3,15)
tether_length = 500
settings = StaticSettings(; rho=1.225, g_earth=[0, 0, -9.806], cd_tether=0.9, d_tether=4,
                    rho_tether=0.85, c_spring=500000)
res!(res, state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings)
```
"""
function res!(res, state_vec, param)
    kite_pos, kite_vel, wind_vel, tether_length, settings, buffers, segments, return_result = param
    par = (kite_pos=kite_pos, kite_vel=kite_vel, wind_vel=wind_vel,
           tether_length=tether_length, settings=settings, segments=segments)
    pj = return_result ? buffers[3] : nothing
    r, T0, p0 = tether_shape(state_vec[1], state_vec[2], state_vec[3], par, pj)
    res .= r
    return_result || return nothing
    return res, MVector(T0), pj, MVector(p0)
end


"""
    get_initial_conditions(filename)

Loads the initialization data for the basic examples and tests. The two
angles in `stateVec` are stored in the MATLAB reference convention and are
converted to this package's elevation/wind-frame-azimuth convention via
[`matlab_to_wind`](@ref) before being returned.

# Arguments
- filename: the filename of the mat file to read

# Returns
- state_vec::MVector{3, Float64} state vector (beta [rad], phi [rad], Tn [N])
  tether orientation and tension at ground station
- kite_pos::MVector{3, Float64} kite position vector in wind reference frame
- kite_vel::MVector{3, Float64} kite velocity vector in wind reference frame
- wind_vel::MMatrix{3, n_nodes(segments), Float64} wind velocity vector in wind reference frame, one column per node
- tether_length: Float64 tether length
- settings::StaticSettings struct containing environmental and tether parameters: see [`StaticSettings`](@ref)
"""
function get_initial_conditions(filename)
    vars = matread(filename)
    sv = vec(get(vars,"stateVec", 0))
    β, φ = matlab_to_wind(sv[1], sv[2])
    state_vec = MVector{3}(β, φ, sv[3])
    kite_pos = MVector{3}(vec(get(vars,"kitePos", 0)))
    kite_vel = MVector{3}(vec(get(vars,"kiteVel", 0)))
    wind_vel = get(vars,"windVel", 0)
    tether_length = get(vars,"tetherLength", 0)

    ENVMT = get(vars,"ENVMT", 0) 
    rho_air = get(ENVMT, "rhos", 0) 
    g_earth = [0; 0; -abs(get(ENVMT, "g", 0))]      # in this way g_earth is a vector [0; 0; -9.81]

    T = get(vars,"T", 0)
    cd_tether = get(T, "CD_tether", 0) 
    d_tether = get(T, "d_tether", 0)*1000           # tether diameter                  [mm]
    E = get(T, "E", 0) 
    A = get(T, "A", 0)
    c_spring = E*A 
    # `rho_t` in the .mat files is the mass per unit length [kg/m], while `StaticSettings` wants
    # a density [kg/m^3] - the model multiplies by the cross section itself. For these
    # fixtures the quotient is 970.0 kg/m^3, the density of Dyneema, which is what makes
    # the interpretation unambiguous. Passing `rho_t` through unconverted made the tether
    # 1/A = 1442 times too light and was the cause of the ~2846 N gap against `T0` that
    # test/test_qsm.jl used to document as unexplained.
    rho_tether = get(T, "rho_t", 0) / A

    # the .mat files store one wind column per node, so there is one segment more than that
    settings = StaticSettings(; rho=rho_air, g_earth, cd_tether, d_tether, rho_tether, c_spring,
                        segments=size(wind_vel, 2) + 1)

    return state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings
end

"""
    init_quasisteady(kite_pos, tether_length; kite_vel = nothing, segments = nothing, wind_vel = nothing, settings = nothing)

Initialize the quasi-steady tether model providing an initial guess for the state vector based on the numerical solution of the catenary equation

# Arguments
- kite_pos::MVector{3, Float64} kite position vector in wind reference frame
- tether_length: Float64 tether length
- kite_vel::MVector{3, Float64} kite velocity vector in wind reference frame
- segments::Int number of tether segments; `segments - 1` nodes are stored
- wind_vel::MMatrix{3, n_nodes(segments), Float64} wind velocity vector in wind reference frame, one column per node
- settings::StaticSettings struct containing environmental and tether parameters: see [`StaticSettings`](@ref)

# Returns
- state_vec::MVector{3, Float64} state vector (beta [rad], phi [rad], Tn [N])  
  tether orientation and tension at ground station
"""
function init_quasisteady(kite_pos, tether_length; kite_vel = nothing, segments = nothing, wind_vel = nothing, settings = nothing)
    # Some basic checks
    @assert isa(kite_pos, MVector{3}) || error("kite_pos must be a MVector of size (3,1)")
    if isnothing(kite_vel) 
        kite_vel = MVector{3}([0.0, 0.0, 0.0])
    end

    # `wind_vel` has one column per node, so one less than `segments`, see `n_nodes`
    if isnothing(segments) && isnothing(wind_vel)
        segments = 8
        wind_vel = zeros(3, n_nodes(segments))
    elseif isnothing(segments) && !isnothing(wind_vel)
        @assert size(wind_vel)[1] == 3 || error("wind_vel should have 3 rows!")
        segments = size(wind_vel)[2] + 1
    elseif !isnothing(segments) && isnothing(wind_vel)
        wind_vel = zeros(3, n_nodes(segments))
    elseif !isnothing(segments) || !isnothing(wind_vel)
        @assert size(wind_vel)[1] == 3 || error("wind_vel should have 3 rows!")
        @assert size(wind_vel)[2] == n_nodes(segments) || error("wind_vel should have n_nodes(segments) columns!")
    end

    if isnothing(settings)
        settings = StaticSettings()
    else
        @assert typeof(settings) == StaticSettings || error("settings should be of type StaticSettings!")
    end

    kite_dist = norm(kite_pos)

    # azimuth angle calculation
    phi_init = atan(kite_pos[2], kite_pos[1])        
    
    function solve_catenary(kite_pos, tether_length, nodes)
        hvec = kite_pos[1:2]    
        h = norm(hvec)
        v = kite_pos[3]
    
        # Function for nonlinear solver
        function f!(res, coeff, param)    
            tether_length, v, h = param  
            res[] = sqrt(tether_length^2 - v^2) - (2 * sinh(coeff[] * h / 2) / coeff[])       
        end
    
        u0 = [0.1]  # Initial guess as a scalar in an array
    
        # Define and solve nonlinear problem
        param = (tether_length, v, h)
        prob = NonlinearProblem(f!, u0, param)
        coeff = solve(prob, NewtonRaphson(autodiff=AutoFiniteDiff()), show_trace=Val(false)) 
        coeff_val = coeff[]  # Extract scalar value
        
        # Adjust catenary solution to specific case
        X = LinRange(0, h, nodes)
        angle1 = atan(hvec[1], hvec[2])
        XY = [sin(angle1) * X'; cos(angle1) * X']
        x_min = -(1 / 2) * (log((tether_length + v) / (tether_length - v)) / coeff_val - h)
        bias = -cosh(-x_min * coeff_val) / coeff_val    
        # Compute z-coordinates of catenary
        z_catenary = cosh.((X .- x_min) .* coeff_val) ./ coeff_val .+ bias
        x_catenary = XY[1, :]
        y_catenary = XY[2, :]    
        return x_catenary, y_catenary, z_catenary, coeff_val
    end

    # Solve the catenary equation
    x_catenary, y_catenary, z_catenary, coeff = solve_catenary(kite_pos, tether_length, n_nodes(segments))
    # Calculate the elevation angle
    beta_init = atan(z_catenary[2], sqrt(x_catenary[2]^2 + y_catenary[2]^2))

    # Initial tension, from the catenary itself: its parameter 1/coeff is H/w, the
    # horizontal tension over the weight per unit length, so the shape that was just
    # fitted to the tether length already carries a tension estimate. The tension of a
    # sagging tether is set by its own weight, not by how stiff it is, which is why the
    # previous guess of a fixed fraction of `c_spring` could be orders of magnitude off -
    # and every decade of error costs `simulate_tether` iterations. `w * tether_length`
    # bounds it from below so that a near vertical tether cannot produce a zero guess.
    w = settings.rho_tether * (π/4 * (settings.d_tether/1000)^2) * abs(settings.g_earth[3])
    tension = sqrt((w / coeff)^2 + (w * tether_length)^2)

    # Assemble state vector
    state_vec = MVector{3}([beta_init, phi_init, tension])
    
    return state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings
end

"""
    get_analytic_catenary(filename)

Loads the analytic catenary curve for the 2D catenary example

# Arguments
- filename: the filename of the mat file to read

# Returns
- x_cat: x coordinates of the catenary curve
- y_cat: x coordinates of the catenary curve
"""
function get_analytic_catenary(filename)
    vars        = matread(filename)
    vars        = get(vars, "analytic_catenary", 0)
    x_cat       = vec(get(vars, "x", 0))
    y_cat       = vec(get(vars, "y", 0))
    return x_cat, y_cat
end

end # module QuasiSteady
