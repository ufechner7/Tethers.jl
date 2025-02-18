using Timers
tic()
using LinearAlgebra
using ModelingToolkit
using OrdinaryDiffEq

using ModelingToolkit: t_nounits as t, D_nounits as D

import Plots
toc()
include("video.jl")

# 🔹 Define struct for simulation settings
G_EARTH::Vector{Float64} = [0.0, -9.81] 
F2::Vector{Float64} = [0.0, 15]
F3::Vector{Float64} = [0.0, 15]
segments::Int64 = 3   
duration::Float64 = 10.0  # Simulation time [s]
save::Bool = false  # Whether to save animation frames
#original positions

POS0 = [
     0.000   -1   1   
     0.000    2   2 
]
VEL0 = zeros(2, segments)


# defining the model, Z component upwards
@parameters k=1.0 m=1.0 l1=sqrt(5) l2=2.0 l3= sqrt(5) damping=0.5
@variables pos(t)[1:2, 1:segments]  = POS0
@variables vel(t)[1:2, 1:segments]  = VEL0
@variables acc(t)[1:2, 1:segments]
@variables segment(t)[1:2, 1:segments]
@variables unit_vector(t)[1:2, 1:segments]
@variables norm1(t)[1:segments]
@variables rel_vel(t)[1:2, 1:segments]
@variables spring_vel(t)[1:segments]
@variables c_spring(t)[1:segments]
@variables spring_force(t)[1:2, 1:segments]
@variables total_force(t)[1:2, 1:segments]


# # 🔹 Define symbolic parameters (Using values from `Settings3`)
# @parameters t k m g[1:2] F2[1:2] F3[1:2] l1 l2 l3

# # 🔹 Define symbolic state variables (positions & velocities)
# # 🔹 Define symbolic state variables (positions & velocities as vectors)
# @variables pos[1:2, 1:3](t)  # Positions (2D) for three points
# @variables vel[1:2, 1:3](t)  # Velocities (2D) for three points
# basic differential equations
eqs1 = vcat(D.(pos) .~ vel,
            D.(vel) .~ acc,
            )
eqs2 = vcat(eqs1...)
eqs2 = vcat(eqs2, pos[:,1] .~ [0.0,0.0])
# # 🔹 Define differential operator
# D = Differential(t)    

# Define which points are connected by each segment
conn = [(1,2), (2,3), (3,1)]  # (P1-P2), (P2-P3), (P3-P1)
# Define rest lengths for each segment
rest_lengths = [l1, l2, l3]  # Each segment has its corresponding rest length

for i in 1:segments  # Loop through each segment
    global eqs2
    local eqs

    eqs = [
        segment[:, i]      ~ pos[:, conn[i][2]] - pos[:, conn[i][1]],  # Compute segment vector
        norm1[i]           ~ norm(segment[:, i]),  # Compute segment length
        unit_vector[:, i]  ~ -segment[:, i] / norm1[i],  # Normalize segment vector
        rel_vel[:, i]      ~ vel[:, conn[i][2]] - vel[:, conn[i][1]],  # Compute relative velocity
        spring_vel[i]      ~ -unit_vector[:, i] ⋅ rel_vel[:, i],  # Project velocity onto spring direction
        
        # Apply Hooke's Law without preventing compression
        c_spring[i]        ~ k,  
        spring_force[:, i] ~ (c_spring[i] * (norm1[i] - rest_lengths[i]) + damping * spring_vel[i]) * unit_vector[:, i]
    ]
    
    eqs2 = vcat(eqs2, reduce(vcat, eqs))
end
 # Loop over only P2 (i=2) and P3 (i=3), since P1 (i=1) is fixed

 for i in 2:3 
    global eqs2
    local eqs = []  # Initialize local equations list for this iteration

    # Compute total force acting on point i (either P2 or P3)
    # Total force is the sum of all forces acting ON this point from connected springs
    incoming_force = sum(
        spring_force[:, j] for j in 1:segments if conn[j][2] == i  # Forces where this point is the end
    )

    outgoing_force = sum(
        spring_force[:, j] for j in 1:segments if conn[j][1] == i  # Forces where this point is the start
    )

    # Apply Newton's Second Law: F = ma, ensuring correct force addition
    if i == 2
        push!(eqs, total_force[:, i] ~ incoming_force - outgoing_force .+ F2)  # Apply F2 as vector addition
    elseif i == 3
        push!(eqs, total_force[:, i] ~ incoming_force - outgoing_force .+ F3)  # Apply F3 as vector addition
    end
    
    push!(eqs, acc[:, i] ~ G_EARTH + total_force[:, i] / m)  # Compute acceleration with gravity

    # Append the equations for this particle to the global system
    eqs2 = vcat(eqs2, reduce(vcat, eqs))
end






# # loop over all tether particles to apply the forces and calculate the accelerations
# for i in 1:(segments)
#     global eqs2; local eqs
#     eqs = []
#     if i == 1
#         push!(eqs, total_force[:, i] ~ spring_force[:, i])
#         push!(eqs, acc[:, i]         ~ zeros(3))
#     else
#         push!(eqs, total_force[:, i] ~ spring_force[:, i-1] - spring_force[:, i])
#         push!(eqs, acc[:, i]         ~ G_EARTH + total_force[:, i] / mass)
#     end
#     eqs2 = vcat(eqs2, reduce(vcat, eqs))
# end




# # loop over all segments to calculate the spring forces
# for i in segments:-1:1
#     global eqs2; local eqs
#     eqs = [segment[:, i]      ~ pos[:, i] - pos[:, i-1],
#            norm1[i]           ~ norm(segment[:, i]),
#            unit_vector[:, i]  ~ -segment[:, i]/norm1[i],
#            rel_vel[:, i]      ~ vel[:, i+1] - vel[:, i],
#            spring_vel[i]      ~ -unit_vector[:, i] ⋅ rel_vel[:, i],
#            c_spring[i]           ~ c_spring0 * (norm1[i] > l_seg),
#            spring_force[:, i] ~ (c_spring[i] * (norm1[i] - l_seg) + damping * spring_vel[i]) * unit_vector[:, i]]
#     eqs2 = vcat(eqs2, reduce(vcat, eqs))
# end

# # 🔹 Define Hooke’s Law for spring force
# function spring_force(xa, ya, xb, yb, l)
#     dx = xb - xa
#     dy = yb - ya
#     length = sqrt(dx^2 + dy^2)  # Compute current length
#     force_mag = k * (length - l)  # Hooke’s Law: F = -k (x - L0)
#     fx = force_mag * dx / length
#     fy = force_mag * dy / length
    
#     return [fx, fy]
# end

# # 🔹 Define fixed anchor point P1 at (0,0)
# pos[:,1] => [0, 0],

# # 🔹 Compute forces symbolically
# F_s1 = spring_force(pos[1,1], pos[2,1], pos[1,2], pos[2,2], l1)  # Force from P1 to P2
# F_s2 = spring_force(pos[1,2], pos[2,2], pos[1,3], pos[2,3], l2)  # Force from P2 to P3
# F_s3 = spring_force(pos[1,3], pos[2,3], pos[1,1], pos[2,1], l3)  # Force from P3 to P1

# # 🔹 Net forces on P2 and P3
# F_net2 = Symbolics.Arr([-F_s1[1] + F_s2[1] + F2[1], -F_s1[2] + F_s2[2] + F2[2]])
# F_net3 = Symbolics.Arr([-F_s2[1] + F_s3[1] + F3[1], -F_s2[2] + F_s3[2] + F3[2]])

# eqs = [
#     D.(pos) .~ vel,
#     D.(vel[:,2]) .~ F_net2 / m,
#     D.(vel[:,3]) .~ F_net3 / m
# ]

@named sys = ODESystem(reduce(vcat, Symbolics.scalarize.(eqs2)), t)
simple_sys = structural_simplify(sys)
# 🔹 Parameter values from `Settings3`
params = Dict(
    m => se.m,
    k => se.k,
    g => se.g_earth,
    F2 => se.F2,
    F3 => se.F3,
    l1 => se.l1,
    l2 => se.l2,
    l3 => se.l3
)

# 🔹 Initial conditions (vectors)
u0 = Dict(
    pos[:,2] => [-1.0, 2.0],
    pos[:,3] => [1.0, 2.0],
    vel[:,2] => [0.0, 0.0],
    vel[:,3] => [0.0, 0.0]
)

# 🔹 Time span
tspan = (0.0, se.duration)

# 🔹 Convert symbolic system to an ODEProblem
prob = ODEProblem(complete(sys), u0, tspan, params)

# 🔹 Solve numerically
sol = solve(prob, Tsit5(), saveat=0.05)

# 🔹 Extract the solution
# 🔹 Extract the solution
pos_sol = sol[pos, :, :]
println(pos_sol)
X = []
Y = []

for i in 1:length(sol.t)  # Iterate over solution time steps
    x = [0, sol[pos[1, 2], i], sol[pos[1, 3], i]]  # Extract X-coordinates
    y = [0, sol[pos[2, 2], i], sol[pos[2, 3], i]]  # Extract Y-coordinates
    push!(X, x)
    push!(Y, y)
end

lines, sc = nothing, nothing
xlim = (-5, 5)
ylim = (0, 6)

for i in 1:length(X)
    x = X[i]
    y = Y[i]
    global lines, sc
    println(x)
    println(y)
    lines, sc = plot_kite(x, y, xlim, ylim, lines, sc)
    # plt.tight_layout()
    plt.pause(0.01)
    plt.show(block=false)
    # sleep(0.05)
end

nothing
