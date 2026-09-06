# A re-usable, acausal tether component, derived from the monolithic model of example
# Tether_08.jl. Instead of one function that builds the whole system, the physics is split
# into components that are wired together with `connect`, following
# https://docs.sciml.ai/ModelingToolkit/stable/basics/Composition/
#
# Components
# ----------
# - `Point3D`      : the connector; across variable `pos`, flow variable `force`
# - `Tether`       : `segments` spring-damper segments with drag, two connectors `p1`, `p2`
# - `FixedEnd`     : holds the point it is attached to at a fixed position
# - `FreeEnd`      : a point mass, falling under gravity and the tether forces
#
# Connection rule
# ---------------
# `Tether` has no inertia at its two end points; it only reports the force it exerts there.
# Therefore **every node must have exactly one component that defines its kinematics**,
# i.e. every `Tether` connector must be connected to a `FixedEnd` or a `FreeEnd` (two
# tethers are joined by connecting both of them to the same `FreeEnd`).
#
# This is a submodule of `Tethers`, so that its names do not end up in `Main`: `runtests.jl`
# includes all examples into the same `Main`, and `Tether_06.jl` defines a global
# `mass_per_meter` there, which a top-level `mass_per_meter(se)` here would collide with.
# Use it with `using Tethers.TetherComponents`.
module TetherComponents

using ModelingToolkit, LinearAlgebra, Parameters
using ModelingToolkit: t_nounits as t, D_nounits as D

# Only distinctive names are exported. `mass_per_meter` and `l_spring` are deliberately not:
# `Tether_06.jl` defines a global `mass_per_meter` in `Main`, and a name exported here that
# collides with one of the other examples breaks `runtests.jl`, which includes all of them
# into the same `Main`. Use `TetherComponents.mass_per_meter(se)` to reach them.
export TetherSettings, set_diameter!, m_end
export Point3D, Tether, FixedEnd, FreeEnd

"""
    TetherSettings

Simulation parameters of the composable tether model; same physics and same units as
`Settings3` of [`Tether_08.jl`](https://github.com/ufechner7/Tethers.jl/blob/main/src/Tether_08.jl).

# Fields
- `g_earth::Vector{Float64}`: gravitational acceleration vector [m/s²]
- `v_wind_tether::Vector{Float64}`: wind velocity acting on the tether [m/s]
- `rho`: air density [kg/m³]
- `cd_tether`: drag coefficient of the tether [-]
- `l0`: initial tether length [m]
- `v_ro`: reel-out speed [m/s]
- `d_tether`: tether diameter [mm]
- `rho_tether`: density of the tether material (Dyneema) [kg/m³]
- `c_spring`: unit spring constant [N]
- `rel_compression_stiffness`: relative compression stiffness [-]
- `damping`: unit damping constant [Ns]
- `segments::Int64`: number of tether segments [-]
- `α0`: initial tether angle [rad]
- `duration`: duration of the simulation [s]
- `save::Bool`: if `true`, save PNG files to the `video` folder
"""
@with_kw mutable struct TetherSettings @deftype Float64
    g_earth::Vector{Float64} = [0.0, 0.0, -9.81] # gravitational acceleration     [m/s²]
    v_wind_tether::Vector{Float64} = [2, 0.0, 0.0]
    rho = 1.225
    cd_tether = 0.958
    l0 = 70                                      # initial tether length             [m]
    v_ro = 0.3                                   # reel-out speed                  [m/s]
    d_tether = 4                                 # tether diameter                  [mm]
    rho_tether = 724                             # density of Dyneema            [kg/m³]
    c_spring = 614600                            # unit spring constant              [N]
    rel_compression_stiffness = 0.01             # relative compression stiffness    [-]
    damping = 473                                # unit damping constant            [Ns]
    segments::Int64 = 10                         # number of tether segments         [-]
    α0 = π/10                                    # initial tether angle            [rad]
    duration = 30                                # duration of the simulation        [s]
    save::Bool = false                           # save png files in folder video
end

"""
    set_diameter!(se, d; c_spring_4mm = 614600, damping_4mm = 473)

Set the tether diameter `d` [mm] in the settings `se` and scale the unit spring constant
`se.c_spring` and the unit damping constant `se.damping` with the cross section, relative
to the reference values of a 4 mm tether.
"""
function set_diameter!(se, d; c_spring_4mm = 614600, damping_4mm = 473)
    se.d_tether = d
    se.c_spring = c_spring_4mm * (d/4.0)^2
    se.damping  = damping_4mm  * (d/4.0)^2
end

"""
    mass_per_meter(se)

Mass per meter of the tether [kg/m], derived from `se.d_tether` and `se.rho_tether`.
"""
mass_per_meter(se) = se.rho_tether * π * (se.d_tether/2000.0)^2

"""
    l_spring(se)

Symbolic (time dependent) length of one tether segment [m], taking reel-out into account.
"""
l_spring(se) = (se.l0 + se.v_ro*t) / se.segments

"""
    m_end(se)

Symbolic (time dependent) mass of a half tether particle [kg], i.e. the tether mass that a
[`FreeEnd`](@ref) attached to a [`Tether`](@ref) built with the same settings has to carry.
"""
m_end(se) = 0.5 * mass_per_meter(se) * l_spring(se)

# ---------------------------------------------------------------------------------------
# connector
# ---------------------------------------------------------------------------------------

"""
    Point3D(; name, pos0=zeros(3))

Connector of a point in 3D space, the 3D equivalent of the `Flange` of the
`Mechanical.Translational` components of the ModelingToolkitStandardLibrary.

Across variable (equal for everything connected to the same node):
- `pos(t)[1:3]`: position [m]

Flow variable (sums up to zero over all components connected to the same node):
- `force(t)[1:3]`: the force that the rest of the system exerts on the component
  owning this connector [N]

The velocity of a node is not part of the connector; a component that needs it uses
`D(pos)`, which keeps the connector balanced (3 across, 3 flow variables).
"""
@connector function Point3D(; name, pos0=zeros(3))
    @variables pos(t)[1:3]
    @variables force(t)[1:3], [connect = Flow]
    # `pos0` is only a guess: the position of a node is determined by the component that
    # defines its kinematics, not by the connector
    guesses = [pos => collect(pos0), force => zeros(3)]
    System(Equation[], t, vcat(collect(pos), collect(force)), []; name, guesses)
end

# ---------------------------------------------------------------------------------------
# the tether itself
# ---------------------------------------------------------------------------------------

"""
    Tether(; name, se, POS0=nothing, VEL0=nothing)

A tether of `se.segments` non-linear spring-damper segments with aerodynamic drag and
reel-out, with the two end points exposed as the connectors `p1` and `p2`.

The tether owns the `se.segments-1` inner particles (full particle mass each) and their
equations of motion. It owns **no** inertia at `p1` and `p2`; there it only reports the
force it exerts, so both connectors must be connected to a [`FixedEnd`](@ref) or a
[`FreeEnd`](@ref) that defines the kinematics of that node and carries the two half
particle masses ([`m_end`](@ref)).

`POS0` and `VEL0` are `3 × (se.segments+1)` matrices of initial positions and velocities of
all particles; they default to the straight line between the connector defaults and zero.

The particle positions are available as the array variable `pos(t)[1:3, 1:se.segments+1]`,
the velocities as `vel`, and the segment lengths as `len(t)[1:se.segments]`.
"""
@component function Tether(; name, se, POS0=nothing, VEL0=nothing)
    n = se.segments
    isnothing(POS0) && (POS0 = zeros(3, n+1))
    isnothing(VEL0) && (VEL0 = zeros(3, n+1))
    size(POS0) == (3, n+1) || error("POS0 must be a 3 x $(n+1) matrix")
    size(VEL0) == (3, n+1) || error("VEL0 must be a 3 x $(n+1) matrix")
    n > 1 || error("a tether needs at least 2 segments, it has $n")

    @named p1 = Point3D(pos0=POS0[:, 1])
    @named p2 = Point3D(pos0=POS0[:, end])

    @parameters rel_compression_stiffness = se.rel_compression_stiffness
    @variables begin
        # the states of this component: the inner particles 2..n
        pos_in(t)[1:3, 1:n-1] = POS0[:, 2:n]
        vel_in(t)[1:3, 1:n-1] = VEL0[:, 2:n]
        acc_in(t)[1:3, 1:n-1]
        # all n+1 particles, including the two end points, for convenient plotting
        pos(t)[1:3, 1:n+1]
        vel(t)[1:3, 1:n+1]
        acc(t)[1:3, 1:n+1]
        segment(t)[1:3, 1:n]
        unit_vector(t)[1:3, 1:n]
        l_seg(t), c_spring(t), damping(t), m_tether_particle(t)
        len(t)[1:n]
        rel_vel(t)[1:3, 1:n]
        spring_vel(t)[1:n]
        c_spr(t)[1:n]
        spring_force(t)[1:3, 1:n]
        v_apparent(t)[1:3, 1:n]
        v_app_perp(t)[1:3, 1:n]
        norm_v_app(t)[1:n]
        half_drag_force(t)[1:3, 1:n]
        total_force(t)[1:3, 1:n+1]
    end

    # the kinematics of the two end particles is defined by whatever is connected to them,
    # the inner particles are integrated here
    eqs2 = vcat([pos[j, 1]   ~ p1.pos[j]    for j in 1:3],
                [vel[j, 1]   ~ D(p1.pos[j]) for j in 1:3],
                [pos[j, n+1] ~ p2.pos[j]    for j in 1:3],
                [vel[j, n+1] ~ D(p2.pos[j]) for j in 1:3],
                [D(pos_in[:, i]) ~ vel_in[:, i] for i in 1:n-1]...,
                [D(vel_in[:, i]) ~ acc_in[:, i] for i in 1:n-1]...,
                # the same particles, but indexed like all the other tether quantities
                [pos[:, i+1] ~ pos_in[:, i] for i in 1:n-1]...,
                [vel[:, i+1] ~ vel_in[:, i] for i in 1:n-1]...,
                [acc[:, i+1] ~ acc_in[:, i] for i in 1:n-1]...)

    # loop over all segments to calculate the spring and drag forces
    for i in 1:n
        eqs = [segment[:, i]      ~ pos[:, i+1] - pos[:, i],
               len[i]             ~ norm(segment[:, i]),
               unit_vector[:, i]  ~ -segment[:, i]/len[i],
               rel_vel[:, i]      ~ vel[:, i+1] - vel[:, i],
               spring_vel[i]      ~ -unit_vector[:, i] ⋅ rel_vel[:, i],
               c_spr[i]           ~ c_spring / (1+rel_compression_stiffness)
                                     * (rel_compression_stiffness+(len[i] > l_seg)),
               spring_force[:, i] ~ (c_spr[i] * (len[i] - l_seg)
                                     + damping * spring_vel[i]) * unit_vector[:, i],
               v_apparent[:, i]   ~ se.v_wind_tether .- (vel[:, i] + vel[:, i+1])/2,
               v_app_perp[:, i]   ~ v_apparent[:, i] - (v_apparent[:, i] ⋅ unit_vector[:, i]) .* unit_vector[:, i],
               norm_v_app[i]      ~ norm(v_app_perp[:, i]),
               half_drag_force[:, i] ~ 0.25 * se.rho * se.cd_tether * norm_v_app[i] * (len[i]*se.d_tether/1000.0)
                                        * v_app_perp[:, i]]
        eqs2 = vcat(eqs2, reduce(vcat, eqs))
    end

    # loop over the inner particles to apply the forces and calculate the accelerations
    for i in 2:n
        eqs = [total_force[:, i] ~ spring_force[:, i-1] - spring_force[:, i] +
                                   half_drag_force[:, i-1] + half_drag_force[:, i],
               acc[:, i]         ~ se.g_earth .+ total_force[:, i] / m_tether_particle]
        eqs2 = vcat(eqs2, reduce(vcat, eqs))
    end

    # the end points carry no mass here, so the force the tether exerts on them
    # (total_force) has to be balanced by the force flowing in through the connector
    eqs = [total_force[:, 1]   ~ -spring_force[:, 1] + half_drag_force[:, 1],
           total_force[:, n+1] ~  spring_force[:, n] + half_drag_force[:, n],
           p1.force            ~ -total_force[:, 1],
           p2.force            ~ -total_force[:, n+1],
           acc[:, 1]           ~ zeros(3),   # unused, the end kinematics is external
           acc[:, n+1]         ~ zeros(3)]
    eqs2 = vcat(eqs2, reduce(vcat, eqs))

    # scalar equations
    eqs = [l_seg             ~ l_spring(se),
           c_spring          ~ se.c_spring / l_seg,
           m_tether_particle ~ mass_per_meter(se) * l_seg,
           damping           ~ se.damping  / l_seg]
    eqs2 = vcat(eqs2, reduce(vcat, eqs))

    # Only `pos_in` and `vel_in` are states and get an initial value (see the `@variables`
    # block above); everything else gets a guess, derived from the shape `POS0`. Without
    # these guesses the initialization of the composed system fails with
    # "Cyclic guesses detected in the system".
    SEG0  = POS0[:, 2:end] - POS0[:, 1:end-1]
    LEN0  = [norm(SEG0[:, i]) for i in 1:n]
    UV0   = -SEG0 ./ LEN0'
    l0    = se.l0 / n
    ZEROS = zeros(3, n)
    VW    = repeat(se.v_wind_tether, 1, n)
    guesses = [pos             => POS0,
               vel             => VEL0,
               acc             => repeat(se.g_earth, 1, n+1),
               acc_in          => repeat(se.g_earth, 1, n-1),
               segment         => SEG0,
               unit_vector     => UV0,
               len             => LEN0,
               rel_vel         => ZEROS,
               spring_vel      => zeros(n),
               c_spr           => fill(se.c_spring/l0, n),
               spring_force    => ZEROS,
               v_apparent      => VW,
               v_app_perp      => VW,
               norm_v_app      => fill(norm(se.v_wind_tether), n),
               half_drag_force => ZEROS,
               total_force     => zeros(3, n+1),
               l_seg           => l0,
               c_spring        => se.c_spring / l0,
               damping         => se.damping / l0,
               m_tether_particle => mass_per_meter(se) * l0]

    System(reduce(vcat, Symbolics.scalarize.(eqs2)), t; name, systems=[p1, p2], guesses)
end

# ---------------------------------------------------------------------------------------
# end point components
# ---------------------------------------------------------------------------------------

"""
    FixedEnd(; name, pos0)

Holds the node it is connected to at the fixed position `pos0`, e.g. a winch or a ground
anchor. The force needed to do so is the force flowing through its connector `flange`.
"""
@component function FixedEnd(; name, pos0)
    @named flange = Point3D(pos0=pos0)
    eqs = collect(flange.pos .~ collect(pos0))
    System(eqs, t; name, systems=[flange])
end

"""
    FreeEnd(; name, se, m_extra=0.0, pos0, vel0=zeros(3))

A point mass at the end of (or between) tethers, free to move under gravity and the forces
flowing in through its connector `flange`.

Its mass is `m_extra` plus one half tether particle mass per tether attached, i.e.
`n_tethers * m_end(se)`; pass the number of tethers connected to it as `n_tethers` and use
the same settings `se` as for those tethers. `m_extra` is the payload mass, e.g. of a kite.
"""
@component function FreeEnd(; name, se, m_extra=0.0, n_tethers=1, pos0, vel0=zeros(3))
    @named flange = Point3D(pos0=pos0)
    @variables begin
        pos(t)[1:3] = pos0
        vel(t)[1:3] = vel0
        acc(t)[1:3]
        mass(t)
    end
    eqs = [pos    ~ flange.pos,
           D(pos) ~ vel,
           D(vel) ~ acc,
           mass   ~ m_extra + n_tethers * m_end(se),
           acc    ~ se.g_earth .+ flange.force / mass]
    guesses = [acc  => copy(se.g_earth),
               mass => m_extra + 0.5 * n_tethers * mass_per_meter(se) * se.l0/se.segments]
    System(reduce(vcat, Symbolics.scalarize.(eqs)), t; name, systems=[flange], guesses)
end

end # module TetherComponents
