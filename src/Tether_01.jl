# Example one: Falling mass.
using Timers
tic()
using ModelingToolkit, OrdinaryDiffEq, ControlPlots
using ModelingToolkit: t_nounits as t, D_nounits as D
toc()

G_EARTH::Vector{Float64} = [0.0, 0.0, -9.81]    # gravitational acceleration     [m/s²]

# defining the model
@variables   pos(t)[1:3] = [0.0, 0.0,  0.0]
@variables   vel(t)[1:3] = [0.0, 0.0, 50.0] 
@variables   acc(t)[1:3]

eqs = vcat(D(pos) ~ vel,
           D(vel) ~ acc,
           acc    ~ G_EARTH)

@named sys = System(eqs, t)
simple_sys = mtkcompile(sys)

# running the simulation
duration = 10.0
dt       = 0.02
tol      = 1e-6
ts       = 0:dt:duration

prob = ODEProblem(simple_sys, nothing, (0.0, duration))
toc()
@time sol = solve(prob, Rodas5(), dt=dt, abstol=tol, reltol=tol, saveat=ts)

# plotting the result
X = sol.t
POS_Z = stack(sol[pos], dims=1)[:,3]
VEL_Z = stack(sol[vel], dims=1)[:,3]

p = plot(X, POS_Z, VEL_Z; xlabel="time [s]", ylabels=["pos_z [m]", "vel_z [m/s]"], 
         labels=["pos_z [m]", "vel_z [m/s]"], fig="falling mass")
display(p)
