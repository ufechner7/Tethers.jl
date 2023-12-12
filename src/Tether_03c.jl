# Example two: Falling mass, attached to non-linear spring without compression stiffness,
# initially moving upwards with 4 m/s.
using ModelingToolkit, OrdinaryDiffEq, PyPlot, LinearAlgebra

G_EARTH  = Float64[0.0, 0.0, -9.81]    # gravitational acceleration [m/s²]
L0 = -10.0                             # initial spring length      [m]
V0 = 4.0                               # initial velocity           [m/s]

function model3(G_EARTH, L0, V0)
    # model, Z component upwards
    @parameters mass=1.0 c_spring0=50.0 damping=0.5 l0=L0
    @variables t pos(t)[1:3] = [0.0, 0.0,  L0]
    @variables   vel(t)[1:3] = [0.0, 0.0,  V0] 
    @variables   acc(t)[1:3] = G_EARTH
    @variables unit_vector(t)[1:3]  = [0.0, 0.0, -sign(L0)]
    @variables c_spring(t) = c_spring0
    @variables spring_force(t)[1:3] = [0.0, 0.0, 0.0]
    @variables force(t) = 0.0 norm1(t) = abs(l0) spring_vel(t) = 0.0
    D = Differential(t)

    eqs = vcat(D.(pos)      ~ vel,
            D.(vel)      ~ acc,
            norm1        ~ norm(pos),
            unit_vector  ~ -pos/norm1,         # direction from point mass to origin
            spring_vel   ~ -unit_vector ⋅ vel,
            c_spring     ~ c_spring0 * (norm1 > abs(l0)),
            spring_force ~ (c_spring * (norm1 - abs(l0)) + damping * spring_vel) * unit_vector,
            acc          ~ G_EARTH + spring_force/mass)

    @named sys = ODESystem(eqs, t)
    simple_sys = structural_simplify(sys)
    simple_sys, pos, vel, c_spring
end

function solve3(simple_sys, L0, V0; cb=true)
    local cb_
    duration = 10.0
    dt = 0.02
    tol = 1e-6
    tspan = (0.0, duration)
    ts    = 0:dt:duration
    # initial state
    u0 = Dict(pos=>[0,0,L0], vel=>[0,0,V0])

    prob = ODEProblem(simple_sys, u0, tspan)
    if cb
        function condition(u, t, integrator) # Event when condition(u,t,integrator) == 0
            norm(u[1:3]) - abs(L0)
        end
        function affect!(integrator) end
        cb_ = ContinuousCallback(condition, affect!; interp_points=2)
        solve(prob, Rodas5(), dt=dt, abstol=tol, reltol=tol, saveat=ts, callback = cb_)
        @time sol = solve(prob, Rodas5(), dt=dt, abstol=tol, reltol=tol, saveat=ts, callback = cb_)
    else
        solve(prob, Rodas5(), dt=dt, abstol=tol, reltol=tol, saveat=ts)
        @time sol = solve(prob, Rodas5(), dt=dt, abstol=tol, reltol=tol, saveat=ts)
    end
    sol
end

function plt(sol; title="")
    fig = figure(title)
    X = sol.t
    POS_Z = sol(X, idxs=pos[3])
    VEL_Z = sol(X, idxs=vel[3])

    lns1 = plot(X, POS_Z, color="green", label="pos_z")
    xlabel("time [s]")
    ylabel("pos_z [m]")
    lns2 = plot(X, L0.+0.005 .* sol[c_spring], color="grey", label="c_spring")
    grid(true)
    twinx()
    ylabel("vel_z [m/s]") 
    lns3 = plot(X, VEL_Z, color="red", label="vel_z")
    lns = vcat(lns1, lns2, lns3)
    labs = [l.get_label() for l in lns]
    legend(lns, labs) 
    PyPlot.show(block=false)
    nothing
end

println("Solving the system without callback...")
simple_sys, pos, vel, c_spring = model3(G_EARTH, L0, V0)
sol = solve3(simple_sys, L0, V0; cb=false)
plt(sol; title="Solution without callback")
println("Press any key...")
if ! @isdefined __PC
    readline()
end
println("Solving the system with callback...")
sol = solve3(simple_sys, L0, V0; cb=true)
plt(sol; title="Solution with callback")
println("If you zoom in to the points in time where pos_z crosses -10m")
println("you should see a difference...")
nothing
