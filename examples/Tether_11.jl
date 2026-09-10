# Tutorial example simulating a 3D mass-spring system with a nonlinear spring (1% stiffness
# for l < l_0), n tether segments, tether drag and reel-in and reel-out.
# New feature: a kite flying a circular trajectory on a cone, as in
# `examples/quasisteady/flying_circular.jl`, but resolved with the dynamic
# mass-spring-damper tether model instead of a quasi-steady catenary. The kite's end of the
# tether is driven by a `MovingEnd` (see `src/TetherComponent.jl`), which imposes the
# trajectory's position as a function of time; the tether component gets the endpoint's
# velocity for free by differentiating that expression, so no separate velocity or
# acceleration needs to be derived by hand.
using ModelingToolkit, OrdinaryDiffEqCore, OrdinaryDiffEqBDF, SteadyStateDiffEq, LinearAlgebra, Timers, MakieControlPlots
tic()
using ModelingToolkit: t_nounits as t
using ADTypes: AutoFiniteDiff
using Tethers: display_if_interactive
using Tethers.TetherComponents: TetherSettings, set_diameter!, FixedEnd, MovingEnd, assemble_tether
import GLMakie

"""
    circular_kite_pos(avg_el, cone_ang, gamma_dot, traj_dist, t)

Position at time `t` of a kite flying a circular trajectory on a cone of half angle
`cone_ang` and radius `traj_dist`, whose axis is tilted by the average elevation `avg_el`
around the x-axis, turning at the constant angular velocity `gamma_dot`. This is the same
parameterization as in `examples/quasisteady/flying_circular.jl`, with `gamma = gamma_dot*t`.

Works both for a plain number `t` and for the symbolic time variable of a ModelingToolkit
model, since it only uses generic arithmetic and trigonometric functions.
"""
function circular_kite_pos(avg_el, cone_ang, gamma_dot, traj_dist, t)
    γ = gamma_dot * t
    s, c = sin(cone_ang), cos(cone_ang)
    pos0 = traj_dist .* [s*cos(γ), s*sin(γ), c]
    rot_mat = [1 0 0; 0 cos(avg_el) -sin(avg_el); 0 sin(avg_el) cos(avg_el)]
    rot_mat*pos0
end

"""
    catenary_positions(p1, p2, L, segments)

The `segments+1` node positions of a tether of unstretched length `L` hanging at rest
between `p1` and `p2` under gravity alone, spaced at equal arc length.

This is the static equilibrium that the mass-spring model settles into, which makes it a
far better starting point for the steady state solver than a straight line: with the 5%
slack of `main` the middle of a 500 m span sags by almost 70 m, and letting the tether fall
that far and then swing itself to rest - on nothing but its own weak aerodynamic drag -
takes more integration steps than the solver is allowed to take.

Returns `nothing` when there is no catenary to compute, i.e. when the tether is not longer
than the straight line between its end points or when it hangs vertically; the caller then
falls back to that straight line.
"""
function catenary_positions(p1, p2, L, segments)
    Δ = p2 - p1
    h = norm(Δ[1:2])                # horizontal span
    v = Δ[3]                        # height difference
    (h < 1e-6 || L <= norm(Δ))    && return nothing
    # A catenary `z(x) = A*cosh((x - x0)/A)` spanning `h` with the height difference `v` has
    # the arc length `2A*sinh(h/(2A)) = sqrt(L² - v²)`. Substituting `u = h/(2A)` turns that
    # into `sinh(u)/u = sqrt(L² - v²)/h`, whose left hand side grows monotonically from 1,
    # so bisection solves it without pulling a nonlinear solver into this example.
    r = sqrt(L^2 - v^2) / h
    lo, hi = 0.0, 1.0
    while sinh(hi)/hi < r
        hi *= 2
    end
    for _ in 1:100
        mid = (lo + hi)/2
        sinh(mid)/mid < r ? (lo = mid) : (hi = mid)
    end
    A  = h / (lo + hi)              # = h/(2u)
    x0 = (h - A*log((L + v)/(L - v))) / 2  # abscissa of the lowest point, from z(h) - z(0) = v
    ê  = Δ[1:2] / h                 # horizontal direction from p1 to p2
    s0 = A*sinh(-x0/A)              # arc length coordinate of p1
    z0 = A*cosh(-x0/A)
    POS = zeros(3, segments+1)
    for i in 1:segments+1
        s = (i-1)/segments * L      # equal arc length steps, because all segments are equally long
        x = x0 + A*asinh((s + s0)/A)
        POS[1:2, i] .= p1[1:2] .+ ê .* x
        POS[3, i]    = p1[3] + A*cosh((x - x0)/A) - z0
    end
    POS
end

"""
    calc_initial_state(p1, p2, L, segments)

`(POS0, VEL0)`: the catenary shape ([`catenary_positions`](@ref)), or a straight line
between `p1` and `p2` if there is no catenary to compute, both at zero velocity.
"""
function calc_initial_state(p1, p2, L, segments)
    VEL0 = zeros(3, segments+1)
    POS0 = catenary_positions(p1, p2, L, segments)
    if isnothing(POS0)
        POS0 = zeros(3, segments+1)
        Δ = (p2-p1) / segments
        for i in 1:segments+1
            POS0[:, i] .= p1 + (i-1) * Δ
        end
    end
    POS0, VEL0
end

"""
    steady_state(se, simple_sys)

Solve `simple_sys` for its steady state and return the tether shape, a
`3 × (se.segments+1)` matrix. Only meaningful with both end points fixed and `se.v_ro`
zero.
"""
function steady_state(se, simple_sys)
    prob = SteadyStateProblem(ODEProblem(simple_sys, nothing, (0.0, se.duration)))
    sol = solve(prob, DynamicSS(FBDF(autodiff=AutoFiniteDiff())); dt=1e-3, abstol=1e-6, reltol=1e-4)
    SciMLBase.successful_retcode(sol) ||
        error("Steady state solver failed with return code $(sol.retcode)!")
    sol.original[simple_sys.tether.pos][end]
end

"""
    model(se; avg_el, cone_ang, gamma_dot, traj_dist)

Build the composed tether model for the settings `se`, while its free end (the kite) flies
a circular trajectory on a cone: half angle `cone_ang` around an axis tilted by the average
elevation `avg_el`, at the constant angular velocity `gamma_dot`, radius `traj_dist`.

Internally, this first builds the model with the kite end fixed at its `t = 0` position and
`se.v_ro` set to zero and solves for the steady-state tether shape, then rebuilds the model
with the original settings, that shape as initial condition, and the kite end driven by a
[`MovingEnd`](@ref) along the trajectory.

Returns `(simple_sys, sys)`.
"""
function model(se; avg_el, cone_ang, gamma_dot, traj_dist)
    p1 = zeros(3)
    p2 = circular_kite_pos(avg_el, cone_ang, gamma_dot, traj_dist, 0.0)
    POS0, VEL0 = calc_initial_state(p1, p2, se.l0, se.segments)
    # find the steady state with the kite end fixed; v_ro must be zero, otherwise there is none
    v_ro = se.v_ro
    se.v_ro = 0
    end1 = FixedEnd(; name=:end1, pos0=p1)
    try
        simple_sys, = assemble_tether(se; end1, end2=FixedEnd(; name=:end2, pos0=p2), POS0, VEL0)
        POS0 = steady_state(se, simple_sys)
    finally
        se.v_ro = v_ro  # restore the reel-out speed, also if the steady state solver failed
    end
    # create the real model, with the steady state shape as initial condition and the kite
    # end driven along its trajectory
    pos_expr = circular_kite_pos(avg_el, cone_ang, gamma_dot, traj_dist, t)  # symbolic in `t`
    assemble_tether(se; end1, end2=MovingEnd(; name=:end2, pos0=p2, pos_expr), POS0, VEL0)
end

"""
    simulate(se, simple_sys)

Simulate the tether model `simple_sys` over the duration `se.duration` with the adaptive
`FBDF` solver, starting with a step size of 0.02s and storing the result on a 0.02s grid.

Returns a tuple `(sol, elapsed_time)` with the `ODESolution` and the elapsed time in
seconds.
"""
function simulate(se, simple_sys)
    dt = 0.02
    tol = 1e-4
    tspan = (0.0, se.duration)
    ts = 0:dt:se.duration
    prob = ODEProblem(simple_sys, nothing, tspan)
    toc()
    elapsed_time = @elapsed sol = solve(prob, FBDF(autodiff=AutoFiniteDiff()); dt, abstol=tol, reltol=tol, saveat=ts)
    SciMLBase.successful_retcode(sol) ||
        error("Simulation failed with return code $(sol.retcode)!")
    sol, elapsed_time
end

"""
    main(; avg_el=deg2rad(70), cone_ang=deg2rad(10), gamma_dot=0.05, traj_dist=nothing)

Build and simulate the tether while its free end (the kite) flies one full revolution of a
circular trajectory on a cone, exactly as parameterized in
`examples/quasisteady/flying_circular.jl`: half angle `cone_ang` around an axis tilted by
the average elevation `avg_el`, at the constant angular velocity `gamma_dot`. `traj_dist`,
the radius of that trajectory, defaults to `se.l0 / 1.05` = 500 m, so that the tether is 5%
slack, exactly as in the quasi-steady example.
"""
function main(; avg_el=deg2rad(70), cone_ang=deg2rad(10), gamma_dot=0.05, traj_dist=nothing)
    global sol, pos, total_force, simple_sys, se
    se = TetherSettings()
    se.v_wind_tether = [0.0, 0.0, 0.0]  # the quasi-steady model runs without wind
    se.l0 = 525
    se.segments = 20
    set_diameter!(se, se.d_tether) # adapt spring and damping constants to tether diameter
    se.duration = 2π / gamma_dot   # one full revolution
    se.v_ro = 0                    # constant tether length while flying the circle
    traj_dist = something(traj_dist, se.l0 / 1.05)

    simple_sys, sys = model(se; avg_el, cone_ang, gamma_dot, traj_dist)
    pos, total_force = simple_sys.tether.pos, simple_sys.tether.total_force
    sol, elapsed_time = simulate(se, simple_sys)
    if @isdefined __PC
        return sol, pos, total_force, simple_sys
    end
    println("Elapsed time: $(elapsed_time) s, speed: $(round(se.duration/elapsed_time)) times real-time")
    sol, pos, total_force, simple_sys
end

# the kite flies one full revolution of a circular trajectory, as in flying_circular.jl
sol, pos, total_force, simple_sys = main();

POS = sol[pos]              # one 3×(segments+1) matrix per saved time step
Ft  = sol[total_force]

show_fig(fig, title) = display_if_interactive(() -> display(GLMakie.Screen(; title), fig))

fig1 = GLMakie.Figure()
ax = GLMakie.Axis3(fig1[1, 1]; title="3D view", xlabel="X [m]", ylabel="Y [m]", zlabel="Z [m]", aspect=:data)
l_tether = GLMakie.scatterlines!(ax, POS[1][1, :], POS[1][2, :], POS[1][3, :])
s_origin = GLMakie.scatter!(ax, [0.0], [0.0], [0.0]; markersize=20, marker=:rect, color=:gray)
s_kite   = GLMakie.scatter!(ax, [POS[1][1, end]], [POS[1][2, end]], [POS[1][3, end]]; markersize=12, marker=:diamond, color=:green)
GLMakie.Legend(fig1[1, 2], [l_tether, s_origin, s_kite], ["Tether", "Origin", "Kite"])
show_fig(fig1, "Initial tether shape")

Ft_kite = reduce(hcat, [F[:, end] for F in Ft])
p_force = plot(sol.t, [Ft_kite[1, :]./1000, Ft_kite[2, :]./1000, Ft_kite[3, :]./1000];
               xlabel="time [s]", ylabel="Force [kN]", labels=["F_x", "F_y", "F_z"],
               title="Tether force components at kite during a circular trajectory",
               fig="Tether force at the kite")
display_if_interactive(p_force)

fig3 = GLMakie.Figure()
ax = GLMakie.Axis3(fig3[1, 1]; title="3D view", xlabel="X [m]", ylabel="Y [m]", zlabel="Z [m]", aspect=:data)
s_origin = GLMakie.scatter!(ax, [0.0], [0.0], [0.0]; markersize=20, marker=:rect, color=:gray)
kite_traj = reduce(hcat, [P[:, end] for P in POS])
l_traj = GLMakie.lines!(ax, kite_traj[1, :], kite_traj[2, :], kite_traj[3, :])
l_tethers = nothing
stride = max(1, length(POS) ÷ 20)  # ~20 tether snapshots spread over the full circle
for ii in 1:stride:length(POS)
    global l_tethers = GLMakie.scatterlines!(ax, POS[ii][1, :], POS[ii][2, :], POS[ii][3, :];
                                             marker=:xcross, color=:orange, linestyle=:dot)
end
GLMakie.Legend(fig3[1, 2], [s_origin, l_traj, l_tethers], ["Origin", "Kite trajectory", "Tethers"])
show_fig(fig3, "Tether shapes along the trajectory")

nothing
