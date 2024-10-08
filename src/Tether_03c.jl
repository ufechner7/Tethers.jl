# Example three: Falling mass, attached to non-linear spring without compression stiffness,
# initially moving upwards with 4 m/s. Comparing results with and without callbacks.
using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra, ControlPlots
using ModelingToolkit: t_nounits as t, D_nounits as D

G_EARTH::Vector{Float64} = [0.0, 0.0, -9.81]    # gravitational acceleration     [m/s²]
L0::Float64 = -10.0                             # initial spring length      [m]
V0::Float64 = 4.0                               # initial velocity           [m/s]

function model3(G_EARTH, L0, V0)
    # defining the model, Z component upwards
    @parameters mass=1.0 c_spring0=50.0 damping=0.5 l0=L0
    @variables   pos(t)[1:3] = [0.0, 0.0,  L0]
    @variables   vel(t)[1:3] = [0.0, 0.0,  V0] 
    @variables   acc(t)[1:3]
    @variables unit_vector(t)[1:3]
    @variables c_spring(t)
    @variables spring_force(t)[1:3]
    @variables force(t) norm1(t) spring_vel(t)

    eqs = vcat(D.(pos)   ~ vel,
            D.(vel)      ~ acc,
            norm1        ~ norm(pos),
            unit_vector  ~ -pos/norm1,         # direction from point mass to origin
            spring_vel   ~ -unit_vector ⋅ vel,
            c_spring     ~ c_spring0 * (norm1 > abs(l0)),
            spring_force ~ (c_spring * (norm1 - abs(l0)) + damping * spring_vel) * unit_vector,
            acc          ~ G_EARTH + spring_force/mass)

    @named sys = ODESystem(Symbolics.scalarize.(reduce(vcat, Symbolics.scalarize.(eqs))), t)
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

function plot2(sol; title="")
    X = sol.t
    POS_Z = stack(sol[pos], dims=1)[:,3]
    VEL_Z = stack(sol[vel], dims=1)[:,3]
    p = plot(X, [POS_Z, L0.+0.005 .* sol[c_spring]], VEL_Z; xlabel="time [s]", ylabels=["pos_z [m]", "vel_z [m/s]"], 
        labels=["pos_z [m]", "c_spring", "vel_z [m/s]"], fig=title)
    display(p)
    nothing
end

println("Solving the system without callback...")
simple_sys, pos, vel, c_spring = model3(G_EARTH, L0, V0)
sol = solve3(simple_sys, L0, V0; cb=false)
plot2(sol; title="Solution without callback")
println("Press any key...")
if ! @isdefined __PC
    readline()
end
println("Solving the system with callback...")
sol = solve3(simple_sys, L0, V0; cb=true)
plot2(sol; title="Solution with callback")
println("If you zoom in to the points in time where pos_z crosses -10m")
println("you should see a difference...")
nothing
