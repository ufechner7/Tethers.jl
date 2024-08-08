# Example two: Falling mass, attached to linear spring
using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra
using ModelingToolkit: t_nounits as t, D_nounits as D

G_EARTH::Vector{Float64} = [0.0, 0.0, -9.81]    # gravitational acceleration     [m/s²]
L0::Float64 = -10.0                             # initial spring length      [m]

# defing the model, Z component upwards
@parameters mass=1.0 c_spring=50.0 damping=0.5 l0=L0
@variables   pos(t)[1:3] = [0.0, 0.0,  L0]
@variables   vel(t)[1:3] = [0.0, 0.0,  0.0] 
@variables   acc(t)[1:3] = G_EARTH
@variables unit_vector(t)[1:3] = [0.0, 0.0, -sign(L0)]
@variables spring_force(t)[1:3] = [0.0, 0.0, 0.0]
@variables spring_vel(t) = 0.0

eqs = vcat(D(pos)      ~ vel,
           D(vel)      ~ acc,
           unit_vector  ~ -pos/norm(pos),         # direction from point mass to origin
           spring_vel   ~ -unit_vector ⋅ vel,
           spring_force ~ (c_spring * (norm(pos) - abs(l0)) + damping * spring_vel) * unit_vector,
           acc          ~ G_EARTH + spring_force/mass)

@named sys = ODESystem(Symbolics.scalarize.(reduce(vcat, Symbolics.scalarize.(eqs))), t)
simple_sys = structural_simplify(sys)
nothing