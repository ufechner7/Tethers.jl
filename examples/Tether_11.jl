# Tutorial example simulating a 3D mass-spring system with a nonlinear spring (1% stiffness
# for l < l_0), n tether segments, tether drag and reel-in and reel-out.
# New feature: if the second extremity is fixed, an acceleration can be prescribed for it
# to simulate the imposed motion of a kite. Here, that motion is a circular trajectory on a
# cone, as in `examples/quasisteady/flying_circular.jl`, but resolved with the dynamic
# mass-spring-damper tether model instead of a quasi-steady catenary.
using ModelingToolkit, OrdinaryDiffEqCore, OrdinaryDiffEqBDF, SteadyStateDiffEq, LinearAlgebra, Timers, Parameters, MakieControlPlots
tic()
using ModelingToolkit: t_nounits as t, D_nounits as D
using LaTeXStrings
using ADTypes: AutoFiniteDiff, AutoForwardDiff
using Tethers: display_if_interactive
# `import`, not `using`: menu.jl runs every example into the same `Main`, and GLMakie
# exports `plot` just like MakieControlPlots, so a `using GLMakie` here would make `plot`
# ambiguous in every example run afterwards.
import GLMakie

@with_kw mutable struct Settings3 @deftype Float64
    g_earth::Vector{Float64} = [0.0, 0.0, -9.81] # gravitational acceleration     [m/s²]
    v_wind_tether::Vector{Float64} = [0.1, 0.0, 0.0]
    rho = 1.225
    cd_tether = 0.958
    l0 = 70                                      # initial tether length             [m]
    v_ro = 0.3                                   # reel-out speed                  [m/s]
    d_tether = 4                                 # tether diameter                  [mm]
    rho_tether = 724                             # density of Dyneema            [kg/m³]
    c_spring = 614600                            # unit spring constant              [N]
    rel_compression_stiffness = 0.01             # relative compression stiffness    [-]
    damping = 473                                # unit damping constant            [Ns]
    segments::Int64 = 6                          # number of tether segments         [-]
    α0 = π/10                                    # initial tether angle            [rad]
    avg_el = deg2rad(70)                         # average elevation of the trajectory cone [rad]
    cone_ang = deg2rad(10)                       # half cone angle of the trajectory  [rad]
    gamma_dot = 0.05                             # angular velocity along the trajectory [rad/s]
    duration = 0                                 # duration of the simulation        [s]
    save::Bool = false                           # save png files in folder video
end

function set_tether_diameter!(se, d; c_spring_4mm = 614600, damping_4mm = 473)
    se.d_tether = d
    se.c_spring = c_spring_4mm * (d/4.0)^2
    se.damping = damping_4mm * (d/4.0)^2
end

"""
    circular_kite_state(se, traj_dist, t)

Position, velocity and acceleration at time `t` of a kite flying a circular trajectory on a
cone of half angle `se.cone_ang` and radius `traj_dist`, whose axis is tilted by the average
elevation `se.avg_el` around the x-axis, turning at the constant angular velocity
`se.gamma_dot`. This is the same parameterization as in
`examples/quasisteady/flying_circular.jl`, with `gamma = se.gamma_dot * t`.

Works both for a plain number `t` and for the symbolic time variable of a ModelingToolkit
model, since it only uses generic arithmetic and trigonometric functions.
"""
function circular_kite_state(se, traj_dist, t)
    γ = se.gamma_dot * t
    s, c = sin(se.cone_ang), cos(se.cone_ang)
    pos0 = traj_dist                .* [s*cos(γ), s*sin(γ), c]
    vel0 = traj_dist*se.gamma_dot   .* [-s*sin(γ), s*cos(γ), 0]
    acc0 = traj_dist*se.gamma_dot^2 .* [-s*cos(γ), -s*sin(γ), 0]
    β = se.avg_el
    rot_mat = [1 0 0; 0 cos(β) -sin(β); 0 sin(β) cos(β)]
    rot_mat*pos0, rot_mat*vel0, rot_mat*acc0
end

function calc_initial_state(se; p1, p2)
    # calculate p2 based on se.α0 and se.l0 if not given
    if isnothing(p2)
        z  = cos(se.α0) * se.l0
        y  = sin(se.α0) * se.l0
        p2 = [p1[1], p1[2] - y, p1[3] - z]
        println("p2: ", p2)
    end
    POS0 = zeros(3, se.segments+1)
    VEL0 = zeros(3, se.segments+1)
    # use a linear interpolation between p1 and p2 for the intermediate points
    for i in 1:se.segments+1
        Δ = (p2-p1) / se.segments
        POS0[:, i] .= p1 + (i-1) * Δ
    end
    POS0, VEL0
end

function model(se; p1=[0,0,0], p2=nothing, fix_p1=true, fix_p2=false, acc_p2 = [0,0,0], vel2=zeros(3))
    if ! isnothing(p1)
        @assert isa(p1, AbstractVector) || error("p1 must be a vector")
        @assert (length(p1) == 3)       || error("p1 must have length 3")
    else
        @assert ! fix_p1                || error("if p1 undefined it cannot be fixed")
    end

    if ! isnothing(p2)
        @assert isa(p2, AbstractVector) || error("p2 must be a vector")
        @assert (length(p2) == 3)       || error("p2 must have length 3")
    else
        @assert ! fix_p2                || error("if p2 undefined it cannot be fixed")
    end
    # straight line approximation for the tether
    POS0, VEL0 = calc_initial_state(se; p1, p2)
    # find steady state
    v_ro = se.v_ro      # save the reel-out speed
    se.v_ro = 0         # v_ro must be zero, otherwise finding the steady state is not possible
    local sol1, pos
    try
        # `acc_p2` has to go for exactly the same reason: with an acceleration prescribed
        # on the second end point, that point never stops accelerating, `norm(du)` never
        # drops below the termination tolerance however loose it is, and there is no
        # steady state to find in the first place. The prescribed acceleration is restored
        # for the real model below, together with the reel-out speed. `vel2` goes for the
        # same reason: a nonzero velocity prescribed on a point held at zero acceleration
        # just translates it forever, so the warm-up pass keeps VEL0 fully zero and only the
        # real model below starts with the prescribed initial velocity.
        simple_sys, pos, =  model(se, p1, p2, true, true, POS0, VEL0, zeros(3))
        tspan = (0.0, se.duration)
        prob = ODEProblem(simple_sys, nothing, tspan)
        prob1 = SteadyStateProblem(prob)
        # the tether swings as a whole for a long time before the per-segment dampers bring
        # it to rest, so DynamicSS's tight default termination tolerance (abstol=1e-8,
        # reltol=1e-6) is never quite met; POS0 is only a warm start for the real
        # simulation below, so a looser tolerance here is fine
        sol1 = solve(prob1, DynamicSS(FBDF(autodiff=AutoFiniteDiff())); dt=1e-3, abstol=1e-6, reltol=1e-4)
    finally
        se.v_ro = v_ro  # restore the reel-out speed, also if the steady state solver failed
    end
    SciMLBase.successful_retcode(sol1) ||
        error("Steady state solver failed with return code $(sol1.retcode)!")
    POS0 = sol1[pos]
    VEL0[:, end] .= vel2  # e.g. the tangential velocity of a kite flying a circular trajectory
    # create the real model
    model(se, p1, p2, fix_p1, fix_p2, POS0, VEL0, acc_p2)
end
function model(se, p1, p2, fix_p1, fix_p2, POS0, VEL0, acc_p2)
    mass_per_meter = se.rho_tether * π * (se.d_tether/2000.0)^2
    @parameters c_spring0=se.c_spring/(se.l0/se.segments) l_seg=se.l0/se.segments
    @parameters rel_compression_stiffness = se.rel_compression_stiffness
    @variables begin 
        pos(t)[1:3, 1:se.segments+1]  = POS0
        vel(t)[1:3, 1:se.segments+1]  = VEL0
        acc(t)[1:3, 1:se.segments+1]
        segment(t)[1:3, 1:se.segments]
        unit_vector(t)[1:3, 1:se.segments]
        l_spring(t), c_spring(t), damping(t), m_tether_particle(t)
        len(t)[1:se.segments]
        rel_vel(t)[1:3, 1:se.segments]
        spring_vel(t)[1:se.segments]
        c_spr(t)[1:se.segments]
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
               c_spr[i]           ~ c_spring / (1+rel_compression_stiffness) 
                                     * (rel_compression_stiffness+(len[i] > l_spring)),
               spring_force[:, i] ~ (c_spr[i] * (len[i] - l_spring) 
                                     + damping * spring_vel[i]) * unit_vector[:, i],
               v_apparent[:, i]   ~ se.v_wind_tether .- (vel[:, i] + vel[:, i+1])/2,
               v_app_perp[:, i]   ~ v_apparent[:, i] - (v_apparent[:, i] ⋅ unit_vector[:, i]) .* unit_vector[:, i],
               norm_v_app[i]      ~ norm(v_app_perp[:, i]),
               half_drag_force[:, i] ~ 0.25 * se.rho * se.cd_tether * norm_v_app[i] * (len[i]*se.d_tether/1000.0)
                                        * v_app_perp[:, i]]
        eqs2 = vcat(eqs2, reduce(vcat, eqs))
    end
    # loop over all tether particles to apply the forces and calculate the accelerations
    for i in 1:(se.segments+1)
        eqs = []
        if i == se.segments+1
            push!(eqs, total_force[:, i] ~ spring_force[:, i-1] + half_drag_force[:, i-1])
            if isnothing(p2) || ! fix_p2
                push!(eqs, acc[:, i]         ~ se.g_earth .+ total_force[:, i] / (0.5 * m_tether_particle))
            else
                push!(eqs, acc[:, i]         ~ acc_p2)
            end
        elseif i == 1
            push!(eqs, total_force[:, i] ~ spring_force[:, i] + half_drag_force[:, i])
            if isnothing(p1) || ! fix_p1
                push!(eqs, acc[:, i]     ~ se.g_earth .+ total_force[:, i] / (0.5 * m_tether_particle))
            else
                push!(eqs, acc[:, i]     ~ zeros(3))                
            end
        else
            push!(eqs, total_force[:, i] ~ spring_force[:, i-1] - spring_force[:, i] 
                                           + half_drag_force[:, i-1] + half_drag_force[:, i])
            push!(eqs, acc[:, i]         ~ se.g_earth .+ total_force[:, i] / m_tether_particle)
        end
        eqs2 = vcat(eqs2, reduce(vcat, eqs))
    end
    # scalar equations
    eqs = [l_spring          ~ (se.l0 + se.v_ro*t) / se.segments,
           c_spring          ~ se.c_spring / l_spring,
           m_tether_particle ~ mass_per_meter * l_spring,
           damping           ~ se.damping  / l_spring]
    eqs2 = vcat(eqs2, reduce(vcat, eqs))  
        
    @named sys = System(reduce(vcat, Symbolics.scalarize.(eqs2)), t)
    simple_sys = mtkcompile(sys)
    simple_sys, pos, vel, len, c_spr, total_force
end

function simulate(se, simple_sys)
    dt = 0.02
    tol = 1e-6
    tspan = (0.0, se.duration)
    ts    = 0:dt:se.duration
    prob = ODEProblem(simple_sys, nothing, tspan)
    toc()
    elapsed_time = @elapsed sol = solve(prob, FBDF(autodiff=AutoForwardDiff()); dt, abstol=tol, reltol=tol, saveat=ts)
    elapsed_time = @elapsed sol = solve(prob, FBDF(autodiff=AutoForwardDiff()); dt, abstol=tol, reltol=tol, saveat=ts)
    sol, elapsed_time
end

"""
    main(; avg_el=deg2rad(70), cone_ang=deg2rad(10), gamma_dot=0.05, traj_dist=nothing)

Build and simulate the tether while its free end (the kite) flies one full revolution of a
circular trajectory on a cone, exactly as parameterized in
`examples/quasisteady/flying_circular.jl`: half angle `cone_ang` around an axis tilted by
the average elevation `avg_el`, at the constant angular velocity `gamma_dot`. `traj_dist`,
the radius of that trajectory, defaults to `0.95 * se.l0` so the (slightly slack) tether can
actually reach it.
"""
function main(; avg_el=deg2rad(70), cone_ang=deg2rad(10), gamma_dot=0.05, traj_dist=nothing)
    global sol, pos, vel, total_force, simple_sys, se
    se = Settings3()
    se.avg_el, se.cone_ang, se.gamma_dot = avg_el, cone_ang, gamma_dot
    set_tether_diameter!(se, se.d_tether) # adapt spring and damping constants to tether diameter
    se.duration = 2π / se.gamma_dot       # one full revolution
    se.v_ro = 0                           # constant tether length while flying the circle
    traj_dist = something(traj_dist, 0.95 * se.l0)

    p1 = [0.0, 0.0, 0.0]
    p2, vel2, = circular_kite_state(se, traj_dist, 0.0)
    acc_p2 = circular_kite_state(se, traj_dist, t)[3]  # symbolic in `t`, re-evaluated every step

    simple_sys, pos, vel, len, c_spr, total_force = model(se; p1, p2, fix_p1=true, fix_p2=true, acc_p2, vel2)
    sol, elapsed_time = simulate(se, simple_sys)
    if @isdefined __PC
        return sol, pos, vel, total_force, simple_sys
    end
    println("Elapsed time: $(elapsed_time) s, speed: $(round(se.duration/elapsed_time)) times real-time")
    sol, pos, vel, total_force, simple_sys
end

# the kite flies one full revolution of a circular trajectory, as in flying_circular.jl
sol, pos, vel, total_force, simple_sys = main();

POS = sol[pos]              # one 3×(segments+1) matrix per saved time step
Ft  = sol[total_force]
gamma_vec = se.gamma_dot .* sol.t

# `display(fig)` re-uses the one GLMakie window, so each figure would replace the previous
# one as soon as it is shown; a fresh `Screen` gives every figure a window of its own, and
# `title` names it, because every window is called "Makie" otherwise. The `Screen` is created
# inside the closure so that nothing opens a window on CI.
show_fig(fig, title) = display_if_interactive(() -> display(GLMakie.Screen(; title), fig))

fig1 = GLMakie.Figure()
ax = GLMakie.Axis3(fig1[1, 1]; title="3D view", xlabel="X [m]", ylabel="Y [m]", zlabel="Z [m]", aspect=:data)
l_tether = GLMakie.scatterlines!(ax, POS[1][1, :], POS[1][2, :], POS[1][3, :])
s_origin = GLMakie.scatter!(ax, [0.0], [0.0], [0.0]; markersize=20, marker=:rect, color=:gray)
s_kite   = GLMakie.scatter!(ax, [POS[1][1, end]], [POS[1][2, end]], [POS[1][3, end]]; markersize=12, marker=:diamond, color=:green)
GLMakie.Legend(fig1[1, 2], [l_tether, s_origin, s_kite], ["Tether", "Origin", "Kite"])
show_fig(fig1, "Initial tether shape")

fig2 = GLMakie.Figure()
ax = GLMakie.Axis(fig2[1, 1]; title="Tether force components at kite during a circular trajectory",
                  xlabel=L"\gamma [rad]", ylabel="Force [kN]")
Ft_kite = reduce(hcat, [F[:, end] for F in Ft])
lx = GLMakie.lines!(ax, gamma_vec, Ft_kite[1, :]./1000)
ly = GLMakie.lines!(ax, gamma_vec, Ft_kite[2, :]./1000)
lz = GLMakie.lines!(ax, gamma_vec, Ft_kite[3, :]./1000)
GLMakie.Legend(fig2[1, 2], [lx, ly, lz], [L"F_x", L"F_y", L"F_z"])
show_fig(fig2, "Tether force at the kite")

fig3 = GLMakie.Figure()
ax = GLMakie.Axis3(fig3[1, 1]; title="3D view", xlabel="X [m]", ylabel="Y [m]", zlabel="Z [m]", aspect=:data)
s_origin = GLMakie.scatter!(ax, [0.0], [0.0], [0.0]; markersize=20, marker=:rect, color=:gray)
kite_traj = reduce(hcat, [P[:, end] for P in POS])
# `lines!`, not `scatter!`: the trajectory has one point per time step, and that many
# overlapping opaque 3D markers lose against each other in the depth test, which drops whole
# stretches of the circle and makes it look like the kite only flies half of it
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
