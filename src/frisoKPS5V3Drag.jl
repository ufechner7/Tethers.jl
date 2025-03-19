#5point model, a point vertical point force at P4 P3 P2 and P5
using Timers
tic()
using LinearAlgebra, ModelingToolkit, OrdinaryDiffEq, ControlPlots 
using ModelingToolkit: t_nounits as t, D_nounits as D
toc()
include("videoKPS5.jl")
# 🔹 Define struct for simulation settings 
# 3D: [x,y,z] , x is the heading , z is up and y perpendicular to both
G_EARTH::Vector{Float64} = [0.0, 0.0, -9.81]
v_wind_tether::Vector{Float64} = [0.0, 100, 0.0]
F1::Vector{Float64} = [0.0, 0.0,  5]
F2::Vector{Float64} = [0.0, 0.0,  22]
F3::Vector{Float64} = [0.0, 0.0, 20]
F4::Vector{Float64} = [0.0, 0,  10]
F5::Vector{Float64} = [0.0,  0, 10]
segments::Int64 = 10   
points::Int64 = 6  
duration::Float64 = 10.0  # Simulation time [s]
save::Bool = false  # Whether to save animation frames
#original positions
POS0 = [
    0.000   1   -1   0    0   0.000
    0.000   0    0   2   -2   0.000
    10.000  13  13   12   12  0.000
 ]
 VEL0 = zeros(3, points)
# defining the model, Z component upwards	
@parameters K1=400 K2=50 K3=200 m=1.0 l1=sqrt(10) l2=2.0 l3= sqrt(10) l4= sqrt(8) l5= sqrt(6) l6= sqrt(6) l7= sqrt(8) l8= sqrt(6) l9= sqrt(6) l10=10 damping=0.9
@parameters rho=1.225 cd_tether=1.28 d_tether=0.01 
@variables pos(t)[1:3, 1:points]  = POS0
@variables vel(t)[1:3, 1:points]  = VEL0
@variables acc(t)[1:3, 1:points]
@variables segment(t)[1:3, 1:segments]
@variables unit_vector(t)[1:3, 1:segments]
@variables norm1(t)[1:segments]
@variables rel_vel(t)[1:3, 1:segments]
@variables spring_vel(t)[1:segments]
@variables c_spring(t)[1:segments]
@variables spring_force(t)[1:3, 1:segments]
@variables v_apparent(t)[1:3, 1:segments]
@variables v_app_perp(t)[1:3, 1:segments]
@variables norm_v_app(t)[1:segments]
@variables half_drag_force(t)[1:3, 1:segments]
@variables total_force(t)[1:3, 1:points]

# basic differential equations
eqs1 = vcat(D.(pos) .~ vel,
          D.(vel) .~ acc,
          )
eqs2 = vcat(eqs1...)
# Fix P6 at (0,0,-10)
eqs2 = vcat(eqs2, pos[:,6] .~ [0.0, 0.0, 0.0])  

# Define which points connect by each segment, their rest lengths and K1 is the tether, K2 are the bridle lines and K3 is the kite
conn = [(1,2), (2,3), (3,1), (1,4), (2,4), (3,4), (1,5), (2,5), (3,5), (1,6)]
rest_lengths = [l1, l2, l3, l4, l5, l6, l7, l8, l9, l10]  
k_segmentss = [K2, K3, K2, K2, K3, K3, K2, K3, K3, K1]

for i in 1:segments  
    global eqs2
    local eqs
    eqs = [
       segment[:, i]      ~ pos[:, conn[i][2]] - pos[:, conn[i][1]],
       norm1[i]           ~ norm(segment[:, i]),
       unit_vector[:, i]  ~ -segment[:, i] / norm1[i], 
       rel_vel[:, i]      ~ vel[:, conn[i][2]] - vel[:, conn[i][1]], 
       spring_vel[i]      ~ -unit_vector[:, i] ⋅ rel_vel[:, i],  
       c_spring[i]        ~ (k_segmentss[i]/(rest_lengths[i])) * (0.25 + 0.75*(norm1[i] > rest_lengths[i])),   #using unit spring constant K, 
       spring_force[:, i] ~ (c_spring[i] * (norm1[i] - rest_lengths[i]) + damping * spring_vel[i]) * unit_vector[:, i],
       v_apparent[:, i]   ~ v_wind_tether .- (vel[:, conn[i][1]] + vel[:, conn[i][2]]) / 2,
       v_app_perp[:, i]   ~ v_apparent[:, i] - (v_apparent[:, i] ⋅ unit_vector[:, i]) .* unit_vector[:, i],
       norm_v_app[i]      ~ norm(v_app_perp[:, i]),
        if !(i in [2, 5, 6, 8, 9])             # Drag only for tether and bridle
            half_drag_force[:, i] ~ 0.25 * rho * cd_tether * norm_v_app[i] * (rest_lengths[i] * d_tether / 1000.0) * v_app_perp[:, i]
        else
            half_drag_force[:, i] ~ zeros(3)
    end
    ]
    eqs2 = vcat(eqs2, reduce(vcat, eqs))
end

# Apply force balance for all points
for i in 1:points  
    global eqs2
    local eqs = []  

    # Compute total force acting on each point
    force = sum([spring_force[:, j] for j in 1:segments if conn[j][2] == i]; init=zeros(3)) -
            sum([spring_force[:, j] for j in 1:segments if conn[j][1] == i]; init=zeros(3)) + 
            sum([half_drag_force[:, j] for j in 1:segments if conn[j][1] == i]; init=zeros(3))+
            sum([half_drag_force[:, j] for j in 1:segments if conn[j][2] == i]; init=zeros(3))

    external_forces = [F1, F2, F3, F4, F5]
    if i <= length(external_forces)
        push!(eqs, total_force[:, i] ~ force + external_forces[i])
    #else        #not sure if this is needed
        #push!(eqs, total_force[:, i] ~ force)
    end
    push!(eqs, acc[:, i] ~ G_EARTH + total_force[:, i] / m)  

    eqs2 = vcat(eqs2, reduce(vcat, eqs))
end

# Solve the system
@named sys = ODESystem(reduce(vcat, Symbolics.scalarize.(eqs2)), t)
simple_sys = structural_simplify(sys) 

# running the simulation
dt = 0.02
tol = 1e-6
tspan = (0.0, duration)
ts    = 0:dt:duration
prob = ODEProblem(simple_sys, nothing, tspan)
elapsed_time = @elapsed sol = solve(prob, Rodas5(); dt, abstol=tol, reltol=tol, saveat=ts) 
# println("Elapsed time: $(elapsed_time) s, speed: $(round(duration/elapsed_time)) times real-time")

#plotting
# Extract solution positions
pos_sol = sol[pos, :]
# Initialize empty lists for all time steps
X, Y, Z = [], [], []
# Extract position data for all points over time
for i in 1:length(pos_sol)
    x, y, z = [], [], []  # Empty lists to store data
    # Extract P1 (First Point)
    push!(x, pos_sol[i][1])  # X coordinate of P1
    push!(y, pos_sol[i][2])  # Y coordinate of P1
    push!(z, pos_sol[i][3])  # Z coordinate of P1
    # Extract all other points (P2 to P6)
    for j in 2:points
       push!(x, pos_sol[i][3*(j-1) + 1])  # X-coordinate of Pj
       push!(y, pos_sol[i][3*(j-1) + 2])  # Y-coordinate of Pj
       push!(z, pos_sol[i][3*(j-1) + 3])  # Z-coordinate of Pj
    end
    # Store time-step data
    push!(X, x)
    push!(Y, y)
    push!(Z, z)
end
# Animation Setup
lines, sc = nothing, nothing
# Set limits for plot, depending on minimum and maximum values of Y, Z
ylim=(minimum(vcat(Y...))-2, maximum(vcat(Y...))+2)
zlim=(minimum(vcat(Z...))-2, maximum(vcat(Z...))+2)
# Iterate through frames and update the plot dynamically
for i in 1:length(Z)
    z = Z[i]
    y = Y[i]
    global lines, sc
    # Update plot with new positions
    lines, sc = plot_kite(y, z, ylim, zlim, lines, sc, conn)
    plt.pause(0.01)  # Pause to create animation effect
    plt.show(block=false)  # Show frame without blocking
end
nothing  # Ensures script runs to completion