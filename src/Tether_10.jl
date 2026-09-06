# Tutorial example simulating a 3D mass-spring system with a nonlinear spring (1% stiffness
# for l < l_0), n tether segments, tether drag and reel-in and reel-out.
# New feature: the physics of example 8 is packaged as a re-usable, acausal component
# (see TetherComponent.jl) with two end points, which is wired to its boundary conditions
# with `connect`, following
# https://docs.sciml.ai/ModelingToolkit/stable/basics/Composition/
#
# `main()`  reproduces example 8: one tether, first end fixed, second end free.
# `main2()` re-uses the same component twice: two tethers of half the length, joined by a
#           point mass, which is only possible because the component is composable.
using ModelingToolkit, OrdinaryDiffEq, SteadyStateDiffEq, LinearAlgebra, Timers, Parameters, MakieControlPlots
using ModelingToolkit: t_nounits as t, D_nounits as D
using ADTypes: AutoFiniteDiff, AutoForwardDiff
using Tethers: display_if_interactive
# the re-usable component, see src/TetherComponent.jl; imported by name, so that including
# this example cannot collide with the globals defined by the other examples
using Tethers.TetherComponents: TetherSettings, set_diameter!, Tether, FixedEnd, FreeEnd

"""
    linear_positions(se; p1, p2)

Compute a linearly interpolated initial position and velocity for each tether particle
between the end points `p1` and `p2`, for the settings `se`.

If `p2` is `nothing`, it is derived from `p1` using `se.α0` and `se.l0`; if `p1` is
`nothing`, it is derived from `p2` the same way. At least one of them must be given.

Returns `(p1, p2, POS0, VEL0)`, with `POS0` and `VEL0` each a `3 × (se.segments+1)` matrix;
`VEL0` is all zeros.
"""
function linear_positions(se; p1, p2)
    isnothing(p1) && isnothing(p2) && error("at least one of p1 and p2 must be defined")
    # calculate the missing end point based on se.α0 and se.l0
    z = cos(se.α0) * se.l0
    y = sin(se.α0) * se.l0
    if isnothing(p2)
        p2 = [p1[1], p1[2] - y, p1[3] - z]
        println("p2: ", p2)
    elseif isnothing(p1)
        p1 = [p2[1], p2[2] + y, p2[3] + z]
        println("p1: ", p1)
    end
    POS0 = zeros(3, se.segments+1)
    VEL0 = zeros(3, se.segments+1)
    # use a linear interpolation between p1 and p2 for the intermediate points
    Δ = (p2 - p1) / se.segments
    for i in 1:se.segments+1
        POS0[:, i] .= p1 + (i-1) * Δ
    end
    p1, p2, POS0, VEL0
end

"""
    build(se; p1, p2, fix_p1, fix_p2, m1=0.0, m2=0.0, POS0, VEL0)

Compose one [`Tether`](@ref) with an end component at each of its two connectors and
return the tuple `(simple_sys, sys, tether)`.

An end point is held in place by a [`FixedEnd`](@ref) if the corresponding `fix_*` flag is
`true`, and by a [`FreeEnd`](@ref) with the payload mass `m1` or `m2` otherwise.
"""
function build(se; p1, p2, fix_p1, fix_p2, m1=0.0, m2=0.0, POS0, VEL0)
    @named tether = Tether(; se, POS0, VEL0)
    end1 = fix_p1 ? FixedEnd(; name=:end1, pos0=p1) :
                    FreeEnd(; name=:end1, se, m_extra=m1, pos0=p1, vel0=VEL0[:, 1])
    end2 = fix_p2 ? FixedEnd(; name=:end2, pos0=p2) :
                    FreeEnd(; name=:end2, se, m_extra=m2, pos0=p2, vel0=VEL0[:, end])
    eqs = [connect(end1.flange, tether.p1),
           connect(tether.p2, end2.flange)]
    @named sys = System(eqs, t; systems=[tether, end1, end2])
    tic()
    simple_sys = mtkcompile(sys)
    toc()
    simple_sys, sys, tether
end

"""
    steady_state(se, simple_sys, p1, p2)

Solve `simple_sys` for its steady state and return the tether shape `POS0`, a
`3 × (se.segments+1)` matrix. Only meaningful if both end points are fixed at `p1` and
`p2` and `se.v_ro` is zero.

Only the inner particles `tether.pos_in` are states of the model; the two end points are
constant and are added back here.
"""
function steady_state(se, simple_sys, p1, p2)
    prob = SteadyStateProblem(ODEProblem(simple_sys, nothing, (0.0, se.duration)))
    sol = solve(prob, DynamicSS(KenCarp4(autodiff=AutoFiniteDiff())))
    SciMLBase.successful_retcode(sol) ||
        error("Steady state solver failed with return code $(sol.retcode)!")
    hcat(p1, sol[simple_sys.tether.pos_in], p2)
end

"""
    model(se; p1=[0,0,0], p2=nothing, fix_p1=true, fix_p2=false, m1=0.0, m2=0.0)

Build the composed tether model for the settings `se`, using a steady-state solve to find
a physically consistent initial tether shape between the end points `p1` and `p2`.

`p1` and `p2` are the two end points of the tether, each a length-3 vector or `nothing`;
an end point that is `nothing` is derived from the other one using `se.α0` and `se.l0`.
`fix_p1` and `fix_p2` control whether the corresponding end point is held fixed; `m1` and
`m2` are the payload masses of the end points that are not fixed.

Internally, this first builds the model with both ends fixed and `se.v_ro` set to zero and
solves for the steady-state tether positions, then rebuilds the model with the original
settings and that shape as initial condition.

Returns `(simple_sys, sys, tether)`.
"""
function model(se; p1=[0,0,0], p2=nothing, fix_p1=true, fix_p2=false, m1=0.0, m2=0.0)
    for (p, fix, name) in ((p1, fix_p1, "p1"), (p2, fix_p2, "p2"))
        if isnothing(p)
            fix && error("if $name undefined it cannot be fixed")
        else
            isa(p, AbstractVector) || error("$name must be a vector")
            length(p) == 3         || error("$name must have length 3")
        end
    end
    p1, p2, POS0, VEL0 = linear_positions(se; p1, p2)
    # find the steady state with both ends fixed; v_ro must be zero, otherwise there is none
    v_ro = se.v_ro
    se.v_ro = 0
    try
        simple_sys, _, _ = build(se; p1, p2, fix_p1=true, fix_p2=true, POS0, VEL0)
        POS0 = steady_state(se, simple_sys, p1, p2)
    finally
        se.v_ro = v_ro  # restore the reel-out speed, also if the steady state solver failed
    end
    # create the real model, with the steady state shape as initial condition
    build(se; p1, p2, fix_p1, fix_p2, m1, m2, POS0, VEL0)
end

"""
    model2(se; p1=[0,0,0], p2=nothing, m_knot=0.0, m2=0.0)

Build a chain of two tethers of half the length `se.l0` each, joined by a [`FreeEnd`](@ref)
with the payload mass `m_knot`. The first end point `p1` is fixed, the last one is free and
carries the payload mass `m2`.

This re-uses the same [`Tether`](@ref) component twice and shows how the knot carries one
half particle mass per tether attached to it (`n_tethers=2`). With `m_knot = 0` and
`m2 = 0` the chain is physically the same system as the single tether of [`model`](@ref):
the segments have the same length, stiffness and damping, and the knot has exactly the
mass of an inner particle, so `main2()` reproduces the result of `main()`.

Returns `(simple_sys, sys, tether1, tether2)`.
"""
function model2(se; p1=[0,0,0], p2=nothing, m_knot=0.0, m2=0.0)
    p1, p2, POS0, VEL0 = linear_positions(se; p1, p2)
    # start from the steady state shape of the equivalent single tether, as in `model`
    v_ro = se.v_ro
    se.v_ro = 0
    try
        ss, _, _ = build(se; p1, p2, fix_p1=true, fix_p2=true, POS0, VEL0)
        POS0 = steady_state(se, ss, p1, p2)
    finally
        se.v_ro = v_ro
    end
    # each of the two tethers spans half of the straight line, and is half as long; both
    # reel out at half the speed, so that the total length is the same as for one tether
    se2 = deepcopy(se)
    se2.l0   = se.l0 / 2
    se2.v_ro = se.v_ro / 2
    n = se.segments
    mid = (n ÷ 2) + 1
    POS_A = POS0[:, 1:mid]
    POS_B = POS0[:, mid:end]
    se2.segments = size(POS_A, 2) - 1
    size(POS_B, 2) - 1 == se2.segments ||
        error("se.segments must be even, so that the two tethers have the same length")

    @named tether1 = Tether(; se=se2, POS0=POS_A)
    @named tether2 = Tether(; se=se2, POS0=POS_B)
    end1  = FixedEnd(; name=:end1, pos0=p1)
    knot  = FreeEnd(; name=:knot, se=se2, m_extra=m_knot, n_tethers=2, pos0=POS0[:, mid])
    end2  = FreeEnd(; name=:end2, se=se2, m_extra=m2, pos0=p2)
    eqs = [connect(end1.flange, tether1.p1),
           connect(tether1.p2, knot.flange, tether2.p1),  # the knot joins both tethers
           connect(tether2.p2, end2.flange)]
    @named sys = System(eqs, t; systems=[tether1, tether2, end1, knot, end2])
    tic()
    simple_sys = mtkcompile(sys)
    toc()
    simple_sys, sys, tether1, tether2
end

"""
    simulate(se, simple_sys)

Simulate the tether model `simple_sys` over the duration `se.duration` with the adaptive
`FBDF` solver, starting with a step size of 0.02s and storing the result on a 0.02s grid.

The solver is run twice; the first call triggers compilation, and only the elapsed time of
the second (timed) call is returned.

Returns a tuple `(sol, elapsed_time)` with the `ODESolution` and the elapsed time in
seconds of the second solve.
"""
function simulate(se, simple_sys)
    dt = 0.02
    tol = 1e-6
    tspan = (0.0, se.duration)
    ts = 0:dt:se.duration
    prob = ODEProblem(simple_sys, nothing, tspan)
    elapsed_time = @elapsed sol = solve(prob, FBDF(autodiff=AutoForwardDiff()); dt, abstol=tol, reltol=tol, saveat=ts)
    elapsed_time = @elapsed sol = solve(prob, FBDF(autodiff=AutoForwardDiff()); dt, abstol=tol, reltol=tol, saveat=ts)
    sol, elapsed_time
end

"""
    play(se, POS; segments=se.segments)

Animate the tether positions `POS` (a vector of `3 × (segments+1)` matrices, one per stored
time step), at half of real-time speed. If `se.save` is `true`, one PNG file per frame is
written to the `video` folder.

Does nothing when running non-interactively or on CI, where there is nothing to look at.
"""
function play(se, POS; segments=se.segments)
    isinteractive() && get(ENV, "CI", "false") == "false" || return nothing
    dt = 0.05
    ylim = (-1.2 * (se.l0 + se.v_ro*se.duration), 0.5)
    xlim = (-se.l0, se.l0)
    se.save && mkpath("video")
    xy = (se.l0/4.2, -7.0)   # text position
    start = time_ns()
    i = 1; j = 0
    for time in 0:dt:se.duration
        # the simulation runs in steps of 20ms, but the plot is updated only every 50ms
        while i < length(POS) && (i-1)*0.02 < time
            i += 1
        end
        display_if_interactive(plot2d, POS[i], time; segments, xlim, ylim, xy)
        se.save && savefig("video/"*"img-"*lpad(j, 4, "0")*".png")
        j += 1
        if time <= dt
            sleep(0.001)
            start = time_ns()
        end
        wait_until(start + 0.5 * time * 1e9)
    end
    se.save && println("Run the script ./bin/export_gif to create the gif file!")
    nothing
end

"""
    main(; p1=[0,0,0], p2=nothing, fix_p1=true, fix_p2=false)

Build, simulate and animate the composed tether model; see [`model`](@ref) for the meaning
of the arguments. Writes the z position and velocity of the last particle to
`output/Tether_10_julia.csv`, so it can be compared with the result of example 8.

Returns `(sol, pos, vel, simple_sys)`.
"""
function main(; p1=[0,0,0], p2=nothing, fix_p1=true, fix_p2=false)
    global sol, pos, vel, simple_sys, sys
    se = TetherSettings()
    set_diameter!(se, se.d_tether) # adapt spring and damping constants to tether diameter
    simple_sys, sys, _ = model(se; p1, p2, fix_p1, fix_p2)
    pos, vel = simple_sys.tether.pos, simple_sys.tether.vel
    sol, elapsed_time = simulate(se, simple_sys)
    # save the z position and velocity of the last particle, for comparison with example 8
    X     = sol.t
    POS_Z = stack(sol[pos], dims=1)[:, 3, se.segments+1]
    VEL_Z = stack(sol[vel], dims=1)[:, 3, se.segments+1]
    mkpath("output")
    open(joinpath("output", "Tether_10_julia.csv"), "w") do io
        println(io, "time,pos_z,vel_z")
        for i in eachindex(X)
            println(io, "$(X[i]),$(POS_Z[i]),$(VEL_Z[i])")
        end
    end
    if @isdefined __PC
        return sol, pos, vel, simple_sys
    end
    play(se, sol[pos])
    println("Elapsed time: $(elapsed_time) s, speed: $(round(se.duration/elapsed_time)) times real-time")
    println("Number of evaluations per step: ", round(sol.stats.nf/(se.duration/0.02), digits=1))
    println("Equations: ", length(equations(sys)))
    println("Equations of simplified system: ", length(equations(simple_sys)))
    sol, pos, vel, simple_sys
end

"""
    main2(; p1=[0,0,0], p2=nothing, m_knot=0.0, m2=0.0)

Build, simulate and animate a chain of two tethers joined by a point mass; see
[`model2`](@ref). Returns `(sol2, simple_sys2)`.
"""
function main2(; p1=[0,0,0], p2=nothing, m_knot=0.0, m2=0.0)
    global sol2, simple_sys2
    se = TetherSettings()
    set_diameter!(se, se.d_tether)
    simple_sys2, sys2, = model2(se; p1, p2, m_knot, m2)
    sol2, elapsed_time = simulate(se, simple_sys2)
    if @isdefined __PC
        return sol2, simple_sys2
    end
    # the two tethers share the knot particle, so drop it from the second one when plotting
    POS = [hcat(a, b[:, 2:end]) for (a, b) in zip(sol2[simple_sys2.tether1.pos],
                                                  sol2[simple_sys2.tether2.pos])]
    play(se, POS)
    println("Elapsed time: $(elapsed_time) s, speed: $(round(se.duration/elapsed_time)) times real-time")
    sol2, simple_sys2
end

# run the simulation with a free (unfixed) second end point, like example 8
main(p2=[-40,0,-47], fix_p2=false);

nothing
