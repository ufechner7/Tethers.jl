# Example three: Falling mass, attached to non-linear spring without compression stiffness,
# initially moving upwards with 4 m/s.
using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra

G_EARTH::Vector{Float64} = [0.0, 0.0, -9.81]    # gravitational acceleration     [m/s²]
L0::Float64 = -10.0                             # initial spring length      [m]
V0::Float64 = 4                                 # initial velocity           [m/s]

# defining the model, Z component upwards
@parameters mass=1.0 c_spring0=50.0 damping=0.5 l0=L0
@variables t pos(t)[1:3] = [0.0, 0.0,  L0]
@variables   vel(t)[1:3] = [0.0, 0.0,  V0] 
@variables   acc(t)[1:3] = G_EARTH
@variables unit_vector(t)[1:3]  = [0.0, 0.0, -sign(L0)]
@variables c_spring(t) = c_spring0
@variables spring_force(t)[1:3] = [0.0, 0.0, 0.0]
@variables force(t) = 0.0 norm1(t) = abs(l0) spring_vel(t) = 0.0
D = Differential(t)

vel = collect(vel)
acc = collect(acc)
pos = collect(pos)

eqs = vcat(D.(pos)      .~ vel,
           D.(vel)      .~ acc,
           norm1        .~ norm(pos),
           unit_vector  .~ -pos/norm1,         # direction from point mass to origin
           spring_vel   .~ -unit_vector ⋅ vel,
           c_spring     .~ c_spring0 * (norm1 > abs(l0)),
           spring_force .~ (c_spring * (norm1 - abs(l0)) + damping * spring_vel) * unit_vector,
           acc          .~ G_EARTH + spring_force/mass)

@named sys = ODESystem(eqs, t)
simple_sys = structural_simplify(sys)

# running the simulation
duration = 10.0
dt = 0.02
tol = 1e-6
tspan = (0.0, duration)
ts    = 0:dt:duration

prob = ODEProblem(simple_sys, nothing, tspan)
@time sol = solve(prob, Rodas5(), dt=dt, abstol=tol, reltol=tol, saveat=ts)

# plotting the result
X = sol.t
POS_Z = stack(sol[pos], dims=1)[:,3]
VEL_Z = stack(sol[vel], dims=1)[:,3]

lns1 = plot(X, POS_Z, color="green", label="pos_z")
xlabel("time [s]")
ylabel("pos_z [m]")
lns2 = plot(X, L0.+0.005 .* sol[c_spring], color="grey", label="c_spring")
PyPlot.grid(true)
legend() 
twinx()
ylabel("vel_z [m/s]") 
lns3 = plot(X, VEL_Z, color="red", label="vel_z")
lns = vcat(lns1, lns2, lns3)
labs = [l.get_label() for l in lns]
legend(lns, labs) 
nothing
