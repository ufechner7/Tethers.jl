# Example three: Falling mass, attached to non-linear spring without compression stiffness,
# initially moving upwards with 4 m/s.
using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra, ControlPlots
using ModelingToolkit: t_nounits as t, D_nounits as D

G_EARTH::Vector{Float64} = [0.0, 0.0, -9.81]    # gravitational acceleration     [m/s²]
L0::Float64 = -10.0                             # initial spring length      [m]
V0::Float64 = 4                                 # initial velocity           [m/s]

# defining the model, Z component upwards
@parameters mass=1.0 c_spring0=50.0 damping=0.5 l0=L0
@variables   pos(t)[1:3] = [0.0, 0.0,  L0]
@variables   vel(t)[1:3] = [0.0, 0.0,  V0] 
@variables   acc(t)[1:3]
@variables unit_vector(t)[1:3]
@variables c_spring(t)
@variables spring_force(t)[1:3]
@variables force(t) norm1(t) spring_vel(t)

eqs = vcat(D(pos)       ~ vel,
           D(vel)       ~ acc,
           norm1        ~ norm(pos),
           unit_vector  ~ -pos / norm1,         # direction from point mass to origin
           spring_vel   ~ -unit_vector ⋅ vel,
           c_spring     ~ c_spring0 * (norm1 > abs(l0)),
           spring_force ~ (c_spring * (norm1 - abs(l0)) + damping * spring_vel) * unit_vector,
           acc          ~ G_EARTH + spring_force/mass)

@named sys = ODESystem((reduce(vcat, Symbolics.scalarize.(eqs))), t)
simple_sys = structural_simplify(sys)

# running the simulation
duration = 10.0
dt       = 0.02
tol      = 1e-6
tspan    = (0.0, duration)
ts       = 0:dt:duration

prob = ODEProblem(simple_sys, nothing, tspan)
@time sol = solve(prob, Rodas5(), dt=dt, abstol=tol, reltol=tol, saveat=ts)

# plotting the result
X = sol.t
POS_Z = stack(sol[pos], dims=1)[:,3]
VEL_Z = stack(sol[vel], dims=1)[:,3]

p = plot(X, [POS_Z, L0.+0.005 .* sol[c_spring]], VEL_Z; xlabel="time [s]", ylabels=["pos_z [m]", "vel_z [m/s]"], 
         labels=["pos_z [m]", "c_spring", "vel_z [m/s]"], fig="falling mass, non-linear spring")
display(p)
