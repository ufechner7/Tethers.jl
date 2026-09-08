# Tutorial example simulating a 3D mass-spring system with n tether segments, tether drag
# and reel-in and reel-out, like Tether_07.jl, but with one difference: the axial force of
# each segment is not a hand-tuned nonlinear spring (1% stiffness for l < l_0) any more, it
# is `analytic_force` from src/analytic_force.jl.
#
# `analytic_force` (derived in docs/segment_force.md) is the tension of a tether that bows
# into a parabola under a transverse load. It is called with its default `segments=Inf`, the
# continuum limit: one segment of this model is a piece of *real* rope hanging between two
# nodes, and that piece does bow under the wind, so the smooth formula is the right local
# law. The `1 - 1/n²` correction describes how a *discretised* tether under-sags relative to
# the smooth curve, which is not what a single physical segment does -- and with `n = 1` it
# would wipe the sag term out altogether, collapsing the formula to a bare `max(0, EA ε)`
# with a discontinuous derivative. The bowing of the tether as a whole is still there on top
# of this, from the node dynamics (many straight segments hinged together), exactly as in
# Tether_07.jl; this only replaces the *local* tension law of a single segment.
#
# Because a bowing segment keeps a little tension from its own sag even when the chord is
# shorter than the unstretched length, the force stays smooth and strictly positive without
# the old 1% `rel_compression_stiffness` hack.
#
# The axial damping is faded out along with the tension by the same formula, via
# `damping_factor` = min(1, analytic_force / |hooke_force|): a slack cable does not damp
# axial motion either, so leaving the damper at full strength kept pumping force into
# segments that carry almost none.
using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra, Timers, Parameters, MakieControlPlots
using ModelingToolkit: t_nounits as t, D_nounits as D
using MakieControlPlots
using ADTypes: AutoForwardDiff
using Tethers: display_if_interactive, analytic_force, damping_factor

@with_kw mutable struct Settings3b @deftype Float64
    g_earth::Vector{Float64} = [0.0, 0.0, -9.81] # gravitational acceleration     [m/s²]
    v_wind_tether::Vector{Float64} = [2, 0.0, 0.0]
    rho = 1.225
    cd_tether = 0.958
    l0 = 50                                      # initial tether length             [m]
    v_ro = 2                                     # reel-out speed                  [m/s]
    d_tether = 4                                 # tether diameter                  [mm]
    rho_tether = 724                             # density of Dyneema            [kg/m³]
    c_spring = 614600                            # unit spring constant              [N]
    damping = 473                                # unit damping constant            [Ns]
    segments::Int64 = 5                          # number of tether segments         [-]
    α0 = π/10                                    # initial tether angle            [rad]
    duration = 10                                # duration of the simulation        [s]
    save::Bool = false                           # save png files in folder video
    damp_mode::Symbol = :ratio                   # EXPERIMENT: :ratio :stiff :none    [-]
    damp_exp = 1.0                               # EXPERIMENT: damp_frac^damp_exp     [-]
    clamp_axial::Bool = false                    # EXPERIMENT: rope cannot push       [-]
    tension_segs = Inf                           # EXPERIMENT: segments= for the tension [-]
end

"""
    seg_tension(v_wind_perp, l_unstretched, l_segment, d_segment, rho, cd_tether, c_spring)

Axial tension [N] of a single segment, via [`analytic_force`](@ref) in its continuum
limit (`segments=Inf`, the default).

Wraps `analytic_force` behind a plain positional signature so it can be `@register_symbolic`d
and used inside the ModelingToolkit equations of [`model`](@ref) below.
"""
function seg_tension(v_wind_perp, l_unstretched, l_segment, d_segment, rho, cd_tether, c_spring,
                     segments)
    se = (; rho, cd_tether, d_tether=d_segment, c_spring)
    analytic_force(se; v_wind_perp, d_segment, l_unstretched, l_segment, segments)
end
@register_symbolic seg_tension(v_wind_perp::Real, l_unstretched::Real, l_segment::Real,
                                d_segment::Real, rho::Real, cd_tether::Real, c_spring::Real,
                                segments::Real)

"""
    seg_damping(v_wind_perp, l_unstretched, l_segment, d_segment, rho, cd_tether, c_spring)

Fraction of the nominal axial damping a single segment still carries, via
[`damping_factor`](@ref) in the same continuum limit; the companion of
[`seg_tension`](@ref) above, and registered for the same reason.
"""
function seg_damping(v_wind_perp, l_unstretched, l_segment, d_segment, rho, cd_tether, c_spring)
    se = (; rho, cd_tether, d_tether=d_segment, c_spring)
    damping_factor(se; v_wind_perp, d_segment, l_unstretched, l_segment)
end
@register_symbolic seg_damping(v_wind_perp::Real, l_unstretched::Real, l_segment::Real,
                                d_segment::Real, rho::Real, cd_tether::Real, c_spring::Real)

"""
EXPERIMENT: fade the damping with the *tangent stiffness* of `analytic_force` instead of
with the force ratio -- classical stiffness-proportional (Kelvin-Voigt) damping, `c ∝ k`.

Differentiating `F³ + a₂F² + a₀ = 0` implicitly, with `a₂ = EA(1 - L/L₀)` and
`a₀ = -w²EA L³/(24 L₀)` (`cn = 1` in the continuum limit), gives the closed form

    dF/dL · L₀/EA = (F² + w²L²/8) / (F (3F + 2a₂))
"""
function seg_damping_stiff(v_wind_perp, l_unstretched, l_segment, d_segment, rho, cd_tether, c_spring)
    se = (; rho, cd_tether, d_tether=d_segment, c_spring)
    f  = analytic_force(se; v_wind_perp, d_segment, l_unstretched, l_segment)
    iszero(f) && return zero(f)
    w  = 0.5 * rho * cd_tether * (d_segment/1000) * v_wind_perp^2
    a2 = c_spring * (1 - l_segment/l_unstretched)
    (f^2 + w^2*l_segment^2/8) / (f * (3f + 2a2))
end
@register_symbolic seg_damping_stiff(v_wind_perp::Real, l_unstretched::Real, l_segment::Real,
                                d_segment::Real, rho::Real, cd_tether::Real, c_spring::Real)

function set_tether_diameter!(se, d; c_spring_4mm = 614600, damping_4mm = 473)
    se.d_tether = d
    se.c_spring = c_spring_4mm * (d/4.0)^2
    se.damping = damping_4mm * (d/4.0)^2
end

function calc_initial_state(se)
    POS0 = zeros(3, se.segments+1)
    VEL0 = zeros(3, se.segments+1)
    for i in 1:se.segments+1
        l0 = -(i-1)*se.l0/se.segments
        v0 = -(i-1)*se.v_ro/se.segments
        POS0[:, i] .= [sin(se.α0) * l0, 0, cos(se.α0) * l0]
        VEL0[:, i] .= [sin(se.α0) * v0, 0, cos(se.α0) * v0]
    end
    POS0, VEL0
end

function model(se)
    POS0, VEL0 = calc_initial_state(se)
    mass_per_meter = se.rho_tether * π * (se.d_tether/2000.0)^2
    @parameters c_spring0=se.c_spring/(se.l0/se.segments) l_seg=se.l0/se.segments
    @parameters damp_exp=se.damp_exp
    @variables begin
        pos(t)[1:3, 1:se.segments+1]  = POS0
        vel(t)[1:3, 1:se.segments+1]  = VEL0
        acc(t)[1:3, 1:se.segments+1]
        segment(t)[1:3, 1:se.segments]
        unit_vector(t)[1:3, 1:se.segments]
        l_spring(t), damping(t), m_tether_particle(t)
        len(t)[1:se.segments]
        rel_vel(t)[1:3, 1:se.segments]
        spring_vel(t)[1:se.segments]
        damp_frac(t)[1:se.segments]
        axial_force(t)[1:se.segments]
        spring_force(t)[1:3, 1:se.segments]
        v_apparent(t)[1:3, 1:se.segments]
        v_app_perp(t)[1:3, 1:se.segments]
        norm_v_app(t)[1:se.segments]
        half_drag_force(t)[1:3, 1:se.segments]
        total_force(t)[1:3, 1:se.segments+1]
    end
    # basic differential equations
    eqs2 = vcat([D(pos[:, i]) ~ vel[:, i] for i in axes(pos, 2)],
                [D(vel[:, i]) ~ acc[:, i] for i in axes(vel, 2)])
    # loop over all segments to calculate the spring and drag forces
    for i in 1:se.segments
        eqs = [segment[:, i]      ~ pos[:, i+1] - pos[:, i],
               len[i]             ~ norm(segment[:, i]),
               unit_vector[:, i]  ~ -segment[:, i]/len[i],
               rel_vel[:, i]      ~ vel[:, i+1] - vel[:, i],
               spring_vel[i]      ~ -unit_vector[:, i] ⋅ rel_vel[:, i],
               v_apparent[:, i]   ~ se.v_wind_tether .- (vel[:, i] + vel[:, i+1])/2,
               v_app_perp[:, i]   ~ v_apparent[:, i] - (v_apparent[:, i] ⋅ unit_vector[:, i]) .* unit_vector[:, i],
               norm_v_app[i]      ~ norm(v_app_perp[:, i]),
               damp_frac[i]       ~ (se.damp_mode === :none ? 1.0 :
                                     ((se.damp_mode === :stiff ? seg_damping_stiff : seg_damping)(
                                         norm_v_app[i], l_spring, len[i], se.d_tether,
                                         se.rho, se.cd_tether, se.c_spring))^damp_exp),
               axial_force[i]     ~ seg_tension(norm_v_app[i], l_spring, len[i], se.d_tether,
                                                se.rho, se.cd_tether, se.c_spring, se.tension_segs)
                                    + damp_frac[i] * damping * spring_vel[i],
               spring_force[:, i] ~ (se.clamp_axial ? max(axial_force[i], 0.0) : axial_force[i])
                                    * unit_vector[:, i],
               half_drag_force[:, i] ~ 0.25 * se.rho * se.cd_tether * norm_v_app[i] * (len[i]*se.d_tether/1000.0)
                                        * v_app_perp[:, i]]
        eqs2 = vcat(eqs2, reduce(vcat, eqs))
    end
    # loop over all tether particles to apply the forces and calculate the accelerations
    for i in 1:(se.segments+1)
        eqs = []
        if i == se.segments+1
            push!(eqs, total_force[:, i] ~ spring_force[:, i-1] + half_drag_force[:, i-1])
            push!(eqs, acc[:, i]         ~ se.g_earth .+ total_force[:, i] / (0.5 * m_tether_particle))
        elseif i == 1
            push!(eqs, total_force[:, i] ~ spring_force[:, i] + half_drag_force[:, i])
            push!(eqs, acc[:, i]         ~ zeros(3))
        else
            push!(eqs, total_force[:, i] ~ spring_force[:, i-1] - spring_force[:, i]
                                           + half_drag_force[:, i-1] + half_drag_force[:, i])
            push!(eqs, acc[:, i]         ~ se.g_earth .+ total_force[:, i] / m_tether_particle)
        end
        eqs2 = vcat(eqs2, reduce(vcat, eqs))
    end
    # scalar equations
    eqs = [l_spring          ~ (se.l0 + se.v_ro*t) / se.segments,
           m_tether_particle ~ mass_per_meter * l_spring,
           damping           ~ se.damping  / l_spring]
    eqs2 = vcat(eqs2, reduce(vcat, eqs))

    @named sys = System(reduce(vcat, Symbolics.scalarize.(eqs2)), t)
    simple_sys = mtkcompile(sys)
    sys, simple_sys, pos, vel
end

function simulate(se, simple_sys)
    dt = 0.02
    tol = 1e-6
    tspan = (0.0, se.duration)
    ts    = 0:dt:se.duration
    prob = ODEProblem(simple_sys, nothing, tspan)
    elapsed_time = @elapsed sol = solve(prob, FBDF(autodiff=AutoForwardDiff()); dt, abstol=tol, reltol=tol, saveat=ts)
    elapsed_time = @elapsed sol = solve(prob, FBDF(autodiff=AutoForwardDiff()); dt, abstol=tol, reltol=tol, saveat=ts)
    sol, elapsed_time
end

function play(se, sol, pos)
    dt = 0.05
    ylim = (-1.2 * (se.l0 + se.v_ro*se.duration), 0.5)
    xlim = (-se.l0/2, se.l0/2)
    mkpath("video")
    z_max = 0.0
    # text position
    xy = (se.l0/4.2, z_max-7)
    start = time_ns()
    i = 1; j = 0
    for time in 0:dt:se.duration
        # while we run the simulation in steps of 20ms, we update the plot only every 150ms
        # therefore we have to skip some steps of the result
        while sol.t[i] < time
            i += 1
        end
        display_if_interactive(plot2d, sol[pos][i], time; segments=se.segments, xlim, ylim, xy)
        if se.save
            savefig("video/"*"img-"*lpad(j, 4, "0")*".png")
        end
        j += 1
        if time <= dt
            sleep(0.001)
            start = time_ns()
        end
        wait_until(start + 0.5 * time * 1e9)
    end
    if se.save
        println("Run the script ./bin/export_gif to create the gif file!")
    end
    nothing
end

function main()
    global sol, pos, vel
    se = Settings3b()
    set_tether_diameter!(se, se.d_tether) # adapt spring and damping constants to tether diameter
    sys, simple_sys, pos, vel = model(se)
    sol, elapsed_time = simulate(se, simple_sys)
    # saving the z position and velocity of the last particle for comparison
    # with the Python implementation
    X     = sol.t
    POS_Z = stack(sol[pos], dims=1)[:, 3, se.segments+1]
    VEL_Z = stack(sol[vel], dims=1)[:, 3, se.segments+1]
    mkpath("output")
    open(joinpath("output", "Tether_07b_julia.csv"), "w") do io
        println(io, "time,pos_z,vel_z")
        for i in eachindex(X)
            println(io, "$(X[i]),$(POS_Z[i]),$(VEL_Z[i])")
        end
    end
    play(se, sol, pos)
    println("Elapsed time: $(elapsed_time) s, speed: $(round(se.duration/elapsed_time)) times real-time")
    println("Number of evaluations per step: ", round(sol.stats.nf/(se.duration/0.02), digits=1))
    sol, pos, vel, sys, simple_sys
end

if (! @isdefined __BENCH__) || __BENCH__ == false
    sol, pos, vel, sys, simple_sys = main()
end
__BENCH__ = false
nothing
