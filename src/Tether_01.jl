# Example one: Falling mass.
using ModelingToolkit, OrdinaryDiffEq, PyPlot

using PyCall
pygui_stop_all()
pygui_start(:qt5)

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
u0 = zeros(6)
u0[6] = 50.0

prob = ODEProblem(simple_sys, u0, tspan)
@time sol = solve(prob, Rodas5(), dt=dt, abstol=tol, reltol=tol, saveat=ts)

X = sol.t
POS_Z = sol(X, idxs=pos[3])
VEL_Z = sol(X, idxs=vel[3])

plot(X, POS_Z, color="green")
xlabel("time [s]")
grid(true)
twinx()
ylabel("vel_z [m/s]") 
plot(X, VEL_Z, color="red") 
PyPlot.show(block=true)
nothing