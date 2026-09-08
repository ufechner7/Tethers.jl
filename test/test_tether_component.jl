# Unit tests of the re-usable tether component (src/TetherComponent.jl), independent of the
# examples. Each test composes one `Tether` with a `FixedEnd` or a `FreeEnd` at each of its
# two end points, and checks the result against what can be calculated analytically:
#
# - steady state          : the segment lengths of a tether hanging from a fixed point
# - drag                  : the drag of a vertical tether in a side wind
# - catenary              : the shape of a tether hanging between two points at z=0
# - compression stiffness : the stiffness of a compressed and of a stretched segment
using Test, LinearAlgebra, ModelingToolkit, OrdinaryDiffEqCore, OrdinaryDiffEqBDF, SteadyStateDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D
using ADTypes: AutoFiniteDiff
using Tethers.TetherComponents: TetherSettings, set_diameter!, Tether, FixedEnd, FreeEnd

# `runtests.jl` includes all test scripts and examples into the same `Main`, and
# `Tether_10.jl` defines the globals `linear_positions`, `build` and `steady_state` there.
# The helpers below therefore have distinctive names, so that including this script cannot
# add a method to, or overwrite, one of those.

"""
    component_positions(se, p1, p2)

Initial particle positions on the straight line from `p1` to `p2`, a
`3 × (se.segments+1)` matrix.
"""
function component_positions(se, p1, p2)
    Δ = (p2 - p1) / se.segments
    reduce(hcat, [p1 + (i-1) * Δ for i in 1:se.segments+1])
end

"""
    build_tether(se, p1, p2; fix_p2, POS0)

Compose one `Tether` with a `FixedEnd` at `p1` and, at `p2`, either a second `FixedEnd`
(`fix_p2=true`) or a `FreeEnd` without payload mass (`fix_p2=false`). `POS0` is the initial
shape of the tether. Returns the simplified system.
"""
function build_tether(se, p1, p2; fix_p2, POS0)
    @named tether = Tether(; se, POS0)
    end1 = FixedEnd(; name=:end1, pos0=p1)
    end2 = fix_p2 ? FixedEnd(; name=:end2, pos0=p2) :
                    FreeEnd(; name=:end2, se, pos0=POS0[:, end])
    eqs = [connect(end1.flange, tether.p1),
           connect(tether.p2, end2.flange)]
    @named sys = System(eqs, t; systems=[tether, end1, end2])
    mtkcompile(sys)
end

"""
    solve_steady_state(se, simple_sys)

Integrate `simple_sys` until it stops changing and return the tether shape, a
`3 × (se.segments+1)` matrix.
"""
function solve_steady_state(se, simple_sys)
    prob = SteadyStateProblem(ODEProblem(simple_sys, nothing, (0.0, se.duration)))
    sol = solve(prob, DynamicSS(FBDF(autodiff=AutoFiniteDiff())))
    SciMLBase.successful_retcode(sol) ||
        error("Steady state solver failed with return code $(sol.retcode)!")
    # `DynamicSS` integrates the model until it stops changing, and `sol.original` is that
    # `ODESolution`. Read the tether shape from it, and not from `sol` itself: a variable
    # that `mtkcompile` turns into an observed one cannot be read from a steady state
    # solution, because the observed function needs the time argument that such a solution
    # does not carry. An `ODESolution` has a time, so reading from it always works.
    sol.original[simple_sys.tether.pos][end]
end

"""
    steady_state_shape(se, p1, p2; fix_p2, POS0, rho_factor=100)

Steady state shape of a tether between `p1` and `p2`, built as in [`build_tether`](@ref)
and solved as in [`solve_steady_state`](@ref), but with the air density of `se` multiplied
by `rho_factor` while it is solved.

What makes the steady state solver slow is not the stretching of the segments, which the
spring dampers kill off quickly, but the swinging of the tether as a whole: the dampers act
along the segments and hardly couple to that, so the tether keeps swinging for hours of
simulated time (69000s and 330000 steps for the catenary below). Aerodynamic drag does damp
the swinging, so a denser atmosphere shortens the transient a lot: with `rho_factor=100`
the two steady state solves of this file get about 13 times faster, 22.0s to 1.6s for the
hanging tether and 15.8s to 1.1s for the catenary.

Without wind this does not change the result: the tether is at rest at its steady state, so
`v_apparent` and with it the drag vanishes there, whatever the air density. Measured, the
segment lengths of the hanging tether move by less than 1e-13m and the shape of the
catenary by less than 2e-6m, which is far below the tolerances used here.
"""
function steady_state_shape(se, p1, p2; fix_p2, POS0, rho_factor=100)
    # with wind the drag does not vanish at the steady state, and a denser atmosphere would
    # push the tether further downwind, i.e. change the very shape that is looked for
    all(iszero, se.v_wind_tether) ||
        error("rho_factor only leaves the steady state unchanged if there is no wind")
    rho = se.rho
    se.rho *= rho_factor
    try
        solve_steady_state(se, build_tether(se, p1, p2; fix_p2, POS0))
    finally
        se.rho = rho   # restore the air density, also if the steady state solver failed
    end
end

"""
    catenary_parameter(span, len)

Parameter `c` of the catenary `z(x) = c * cosh(x/c)` of a chain of the length `len`, hanging
between two points that are `span` apart and at the same height. It is the solution of
`len = 2c * sinh(span/2c)`, found by bisection; `f(c) = 2c * sinh(span/2c) - len` falls
monotonically from `∞` at `c → 0` to `span - len < 0` for large `c`.
"""
function catenary_parameter(span, len)
    len > span || error("the chain must be longer than the span, else it does not sag")
    f(c) = 2c * sinh(span/(2c)) - len
    lo, hi = 1e-6, 1e6
    for _ in 1:100
        mid = (lo + hi) / 2
        f(mid) > 0 ? (lo = mid) : (hi = mid)
    end
    (lo + hi) / 2
end

"""
    initial_values(simple_sys, vars)

Values of the symbolic variables `vars` of `simple_sys` at `t=0`, i.e. at the initial
condition of the model. Solving over the empty time span `(0.0, 0.0)` is the simplest way
to get at the observed variables of the system: it initializes the model and returns the
initial state, without taking a single step. Do not integrate instead; a model that is far
from its equilibrium just fails with `Unstable`, even though its state at `t=0` is fine.
"""
function initial_values(simple_sys, vars)
    prob = ODEProblem(simple_sys, nothing, (0.0, 0.0))
    sol = solve(prob, FBDF())
    SciMLBase.successful_retcode(sol) ||
        error("Solver failed with return code $(sol.retcode)!")
    [sol[var][1] for var in vars]
end

@testset "TetherComponent, steady state" begin
    se = TetherSettings()
    set_diameter!(se, se.d_tether)
    se.v_wind_tether = zeros(3)  # no wind, so that gravity alone shapes the tether
    se.v_ro = 0.0                # a steady state exists only without reel-out

    p1 = [0.0, 0.0, 0.0]         # upper end point, held by the FixedEnd
    p2 = [-40.0, 0.0, -47.0]     # lower end point, less than se.l0 away, so the tether sags
    POS0 = component_positions(se, p1, p2)

    # First find the steady state with both ends fixed; this is a well conditioned problem
    # and gives a shape that is a good starting point for the tether with the free end.
    POS_ff = steady_state_shape(se, p1, p2; fix_p2=true, POS0)
    # Then release the lower end and let it settle under gravity.
    POS = steady_state_shape(se, p1, p2; fix_p2=false, POS0=POS_ff)

    @test size(POS) == (3, se.segments+1)
    @test all(isfinite, POS)

    # the FixedEnd must hold the upper end point at p1
    @test POS[:, 1] ≈ p1 atol=1e-6

    # the FreeEnd must have fallen: it hangs below its initial position and the tether
    # must be taut, so it ends up about se.l0 away from the upper end point
    @test POS[3, end] < POS_ff[3, end]
    @test norm(POS[:, end] - p1) ≈ se.l0 rtol=1e-2

    seg_lengths = [norm(POS[:, i+1] - POS[:, i]) for i in 1:se.segments]
    l_seg = se.l0 / se.segments

    # all segments are stretched: they carry the weight of the tether below them
    @test all(l -> l > l_seg, seg_lengths)

    # Segment 1 is the one attached to the FixedEnd, and it carries the weight of the whole
    # tether hanging below it. Each following segment carries the weight of one particle
    # less, so it is under less tension and stretches less: the steady state segment length
    # must decrease with the segment number.
    @test issorted(seg_lengths, rev=true)
    @test seg_lengths[1] > seg_lengths[end]
end

@testset "TetherComponent, drag" begin
    se = TetherSettings()
    set_diameter!(se, se.d_tether)
    se.v_wind_tether = [8.0, 0.0, 0.0]   # side wind, blowing in x direction
    se.v_ro = 0.0                        # no reel-out, so that l_seg stays constant

    # Both end points are fixed and have the same x and y coordinate, and they are exactly
    # se.l0 apart, so the tether is a straight vertical line of segments of the nominal
    # length. The wind is then exactly perpendicular to every segment, and because the
    # tether is at rest at t=0, the apparent wind is the wind itself. Under these conditions
    # the drag of a segment is the textbook drag of a cylinder in a perpendicular flow.
    p1 = [0.0, 0.0, 0.0]
    p2 = [0.0, 0.0, -se.l0]
    POS0 = component_positions(se, p1, p2)
    simple_sys = build_tether(se, p1, p2; fix_p2=true, POS0)

    len, vel, half_drag = initial_values(simple_sys,
                                         (simple_sys.tether.len, simple_sys.tether.vel,
                                          simple_sys.tether.half_drag_force))
    l_seg = se.l0 / se.segments
    @test len ≈ fill(l_seg, se.segments)   # the tether is straight and unstretched at t=0
    @test all(iszero, vel)                 # and it is at rest, so v_apparent is the wind

    # analytic drag of a cylinder of the length l_seg and the diameter se.d_tether [mm] in
    # a perpendicular flow: F = 0.5 * rho * cd * A * v², with the reference area A = l * d
    v_wind    = norm(se.v_wind_tether)
    area      = l_seg * se.d_tether/1000
    drag_seg  = 0.5 * se.rho * se.cd_tether * area * v_wind^2
    # the model splits the drag of a segment evenly over the two particles it connects, so
    # `half_drag_force` must be half of the drag of the whole segment
    @test all(≈(drag_seg/2), half_drag[1, :])

    # the drag of a perpendicular flow pushes straight downwind, here in x direction
    @test all(iszero, half_drag[2, :])
    @test all(iszero, half_drag[3, :])

    # the drag of the whole tether is that of a cylinder of the length se.l0; summing the
    # half drag forces of all segments twice adds up both halves of each of them
    drag_total = 0.5 * se.rho * se.cd_tether * (se.l0 * se.d_tether/1000) * v_wind^2
    @test 2 * sum(half_drag[1, :]) ≈ drag_total
end

@testset "TetherComponent, catenary" begin
    se = TetherSettings()
    set_diameter!(se, se.d_tether)
    se.v_wind_tether = zeros(3)   # no wind, so that gravity alone shapes the tether
    se.v_ro = 0.0                 # no reel-out, otherwise there is no steady state

    # Both end points are fixed at the same height z=0, and closer to each other than the
    # tether is long, so the tether sags: a chain that hangs under its own weight takes the
    # shape of a catenary, z(x) = c * cosh(x/c).
    span = 50.0
    p1 = [-span/2, 0.0, 0.0]
    p2 = [ span/2, 0.0, 0.0]
    POS0 = component_positions(se, p1, p2)
    POS  = steady_state_shape(se, p1, p2; fix_p2=true, POS0)

    # the FixedEnds must hold the two end points, and the tether must hang in the x-z plane
    @test POS[:, 1]   ≈ p1 atol=1e-6
    @test POS[:, end] ≈ p2 atol=1e-6
    @test maximum(abs.(POS[2, :])) < 1e-8

    # both end points are at the same height, so the shape must be symmetric
    @test POS[3, :] ≈ reverse(POS[3, :]) atol=1e-6

    # The analytic catenary is that of an inextensible chain, so compare it with a tether
    # of the length it actually has: the segments are stretched by their own weight, but by
    # less than 1e-5 of their length, which is why the comparison with a catenary works.
    len = sum(norm(POS[:, i+1] - POS[:, i]) for i in 1:se.segments)
    @test len ≈ se.l0 rtol=1e-3

    c = catenary_parameter(span, len)
    z_catenary(x) = c * (cosh(x/c) - cosh(span/(2c)))  # zero at the two end points
    sag = c * (cosh(span/(2c)) - 1)                    # depth of the lowest point [m]

    # The tether is a chain of se.segments straight segments with lumped masses, so its
    # particles follow the catenary of the continuous chain only approximately; the error
    # is quadratic in the segment length (0.6 % of the sag for 10, 0.15 % for 20 segments).
    @test maximum(abs.(POS[3, :] - z_catenary.(POS[1, :]))) < 0.01 * sag

    # the deepest point of the tether must match the analytically calculated sag
    @test -minimum(POS[3, :]) ≈ sag rtol=0.01
end

@testset "TetherComponent, compression stiffness" begin
    se = TetherSettings()
    set_diameter!(se, se.d_tether)
    se.v_wind_tether = zeros(3)   # no drag, the spring force alone is of interest here
    se.v_ro = 0.0                 # no reel-out, so that l_seg stays constant here

    # A vertical tether whose segments are alternately stretched and compressed by the same
    # amount. The two end points are fixed, and because the stretched and the compressed
    # segments cancel out, they are exactly se.l0 apart.
    l_seg = se.l0 / se.segments
    Δ = 0.1 * l_seg
    lengths = [isodd(i) ? l_seg + Δ : l_seg - Δ for i in 1:se.segments]
    POS0 = zeros(3, se.segments+1)
    POS0[3, :] = vcat(0.0, -cumsum(lengths))
    p1, p2 = POS0[:, 1], POS0[:, end]
    @test p2[3] ≈ -se.l0

    simple_sys = build_tether(se, p1, p2; fix_p2=true, POS0)
    len, c_spr, c_spring, spring_force =
        initial_values(simple_sys, (simple_sys.tether.len, simple_sys.tether.c_spr,
                                    simple_sys.tether.c_spring,
                                    simple_sys.tether.spring_force))
    @test len ≈ lengths            # the segments are stretched/compressed as prescribed

    # A segment that is longer than l_seg is stretched and has the full stiffness, a shorter
    # one is compressed and only has the relative compression stiffness. The model scales
    # both with 1/(1+rel_compression_stiffness), so the compressed segment ends up with
    # c_spring * rcs/(1+rcs) and the stretched one with the full c_spring.
    rcs = se.rel_compression_stiffness
    @test all(≈(c_spring), c_spr[1:2:end])                 # stretched segments
    @test all(≈(c_spring * rcs/(1+rcs)), c_spr[2:2:end])   # compressed segments

    # the point of the exercise: a compressed segment is only about 1% as stiff as a
    # stretched one, so that the tether can go slack but hardly pushes
    stiffness_ratio = c_spr[2] / c_spr[1]
    @test stiffness_ratio ≈ rcs/(1+rcs)
    @test stiffness_ratio ≈ 0.01 rtol=0.02

    # the same ratio must show up in the force: the segments are stretched and compressed by
    # the same Δ, and the tether is at rest, so the damping term is zero and the spring
    # force is proportional to the stiffness
    f_stretched  = abs.(spring_force[3, 1:2:end])
    f_compressed = abs.(spring_force[3, 2:2:end])
    @test all(≈(c_spring * Δ), f_stretched)
    @test all(≈(stiffness_ratio), f_compressed ./ f_stretched)

    # a stretched segment pulls its end points together, a compressed one pushes them apart
    @test all(spring_force[3, 1:2:end] .* spring_force[3, 2:2:end] .< 0)
end
nothing
