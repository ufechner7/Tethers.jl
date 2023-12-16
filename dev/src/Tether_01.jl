# Example one: Falling mass.
using ModelingToolkit, OrdinaryDiffEq

G_EARTH  = Float64[0.0, 0.0, -9.81]    # gravitational acceleration

# model
@variables t pos(t)[1:3]=[0.0, 0.0,  0.0]
@variables   vel(t)[1:3]=[0.0, 0.0, 50.0] 
@variables   acc(t)[1:3]=[0.0, 0.0, -9.81] 
D = Differential(t)

eqs = vcat(D.(pos) ~ vel,
           D.(vel) ~ acc,
           acc    .~ G_EARTH)

@named sys = ODESystem(eqs, t)
simple_sys = structural_simplify(sys)

duration = 10.0
dt = 0.02
tol = 1e-6
tspan = (0.0, duration)
ts    = 0:dt:duration

# initial state
u0 = Dict(vel=>[0, 0, 50.0])

prob = ODEProblem(simple_sys, u0, tspan)
@time sol = solve(prob, Rodas5(), dt=dt, abstol=tol, reltol=tol, saveat=ts)

X = sol.t
POS_Z = stack(sol[pos], dims=1)[:,3]
VEL_Z = stack(sol[vel], dims=1)[:,3]

plot(X, POS_Z, color="green")
xlabel("time [s]")
ylabel("pos_z [m]")
PyPlot.grid(true)
twinx()
ylabel("vel_z [m/s]") 
plot(X, VEL_Z, color="red")
nothing