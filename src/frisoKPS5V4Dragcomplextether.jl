using Timers
tic()
using LinearAlgebra, ModelingToolkit, OrdinaryDiffEq, ControlPlots 
using ModelingToolkit: t_nounits as t, D_nounits as D
toc()
include("videoKPS5.jl")

# 🔹 Define struct for simulation settings 
# 3D: [x,y,z] , x is the heading , z is up and y perpendicular to both
G_EARTH::Vector{Float64} = [0.0, 0.0, -9.81]
v_wind_tether::Vector{Float64} = [0.0, 15 , 0.0]
F1::Vector{Float64} = [0.0, 0.0,  0.0]
F2::Vector{Float64} = [0.0, 0.0,  192]
F3::Vector{Float64} = [0.0, 0.0, 192]
F4::Vector{Float64} = [0.0, 0,  192]
F5::Vector{Float64} = [0.0,  0, 192]
tethersegments::Int64 = 12
segments::Int64 = 9 + tethersegments    
points::Int64 = 5 + tethersegments  
duration::Float64 = 10.0  # Simulation time [s]
save::Bool = false  # Whether to save animation frames

# Original positions
POS0 = [0.000   1    -1    0    0    0.000;
        0.000   0     0    2   -2    0.000;
        10.000  13   13   12   12    0.000]
if tethersegments > 1
    extra_nodes = [POS0[:,6] + (POS0[:,1] - POS0[:,6]) * i / tethersegments for i in 1:(tethersegments-1)]
    POS0 = hcat(POS0, extra_nodes...)
end     
VEL0 = zeros(3, points)

#add diferent masses for kite, total = 6.2 (KCU 8.4) at bridle
 

# defining the model, Z component upwards	
@parameters K1=614600  K2=10000 K3=60000 m_kite=6.2/4 m_bridle=8.4 rho_tether=724 l1=sqrt(10) l2=2.0 l3=sqrt(10) l4=sqrt(8) l5=sqrt(6) l6=sqrt(6) l7=sqrt(8) l8=sqrt(6) l9=sqrt(6) l10=10 damping=0.9
@parameters rho=1.225 cd_tether=0.958 d_tether=0.004
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
# New observable variable to record drag force for each segment:
@variables drag_force(t)[1:3, 1:segments]
@variables total_force(t)[1:3, 1:points]

# basic differential equations for positions and velocities
eqs1 = vcat(D.(pos) .~ vel,
            D.(vel) .~ acc)
eqs2 = vcat(eqs1...)
# Fix P6 at (0,0,0) origin
eqs2 = vcat(eqs2, acc[:,6] .~ [0.0, 0.0, 0.0])

# Define which points connect by each segment, their rest lengths and stiffness constants.
conn = [(1,2), (2,3), (3,1), (1,4), (2,4), (3,4), (1,5), (2,5), (3,5)]
conn = vcat(conn, [(6+i, 6+i+1) for i in 0:(tethersegments-2)]...)      # adding tether segments connections
conn = vcat(conn, [(6+tethersegments-1, 1)])
rest_lengths = [l1, l2, l3, l4, l5, l6, l7, l8, l9]
rest_lengths = vcat(rest_lengths, [l10/tethersegments for _ in 1:tethersegments]...)
k_segments = [K2, K3, K2, K2, K3, K3, K2, K3, K3]
k_segments = vcat(k_segments, [K1 for _ in 1:tethersegments]...)

# Define equations for each segment including drag force
for i in 1:segments  
    global eqs2
    local eqs = [
       segment[:, i]      ~ pos[:, conn[i][2]] - pos[:, conn[i][1]],
       norm1[i]           ~ norm(segment[:, i]),
       unit_vector[:, i]  ~ -segment[:, i] / norm1[i], 
       rel_vel[:, i]      ~ vel[:, conn[i][2]] - vel[:, conn[i][1]], 
       spring_vel[i]      ~ -unit_vector[:, i] ⋅ rel_vel[:, i],  
       c_spring[i]        ~ (k_segments[i]/rest_lengths[i]) * (0.1 + 0.9*(norm1[i] > rest_lengths[i])),
       spring_force[:, i] ~ (c_spring[i]*(norm1[i] - rest_lengths[i]) + damping * spring_vel[i]) * unit_vector[:, i],
       v_apparent[:, i]   ~ v_wind_tether .- (vel[:, conn[i][1]] + vel[:, conn[i][2]]) / 2,
       v_app_perp[:, i]   ~ v_apparent[:, i] - (v_apparent[:, i] ⋅ unit_vector[:, i]) .* unit_vector[:, i],
       norm_v_app[i]      ~ norm(v_app_perp[:, i]),
       half_drag_force[:, i] ~ 0.25 * rho * cd_tether * norm_v_app[i] * (rest_lengths[i]*d_tether) * v_app_perp[:, i]#,
       #drag_force[:, i]   ~ half_drag_force[:, i]  # assign computed drag to drag_force
    ]
    eqs2 = vcat(eqs2, reduce(vcat, eqs))
end
#mass distribution
mass_tether = (d_tether^2 )* pi * rho_tether*l10
mass_tetherpoints = mass_tether/(tethersegments+1)
PointMassesBridleKite = [m_bridle+mass_tetherpoints, m_kite, m_kite, m_kite, m_kite]
PointMasses = vcat(PointMassesBridleKite, [mass_tetherpoints for _ in 1:tethersegments]...)
# Apply force balance for all points
for i in 1:points  
    global eqs2
    local eqs = []  
    # Compute total force acting on each point from spring forces and drag forces.
    force = sum([spring_force[:, j] for j in 1:segments if conn[j][2] == i]; init=zeros(3)) -
            sum([spring_force[:, j] for j in 1:segments if conn[j][1] == i]; init=zeros(3)) +
            sum([half_drag_force[:, j] for j in 1:segments if conn[j][1] == i]; init=zeros(3)) +
            sum([half_drag_force[:, j] for j in 1:segments if conn[j][2] == i]; init=zeros(3))
    ExternalForces = [F1, F2, F3, F4, F5]
    if i <= length(ExternalForces)
        push!(eqs, total_force[:, i] ~ force + ExternalForces[i])
    elseif i != 6        # (optional additional constraint)                  
        push!(eqs, total_force[:, i] ~ force)
    end
    push!(eqs, acc[:, i] ~ G_EARTH + total_force[:, i] / PointMasses[i])
    eqs2 = vcat(eqs2, reduce(vcat, eqs))
end

# Solve the system
@named sys = ODESystem(reduce(vcat, Symbolics.scalarize.(eqs2)), t)
simple_sys = structural_simplify(sys) 

dt = 0.02
tol = 1e-6
tspan = (0.0, duration)
ts    = 0:dt:duration
prob = ODEProblem(simple_sys, nothing, tspan)
elapsed_time = @elapsed sol = solve(prob, Rodas5(); dt, abstol=tol, reltol=tol, saveat=ts) 
println("Elapsed time: $(elapsed_time) s, speed: $(round(duration/elapsed_time)) times real-time")

# # Extract the total force acting on each point over time
# spring = sol[spring_force, :]

#Print the total force for each point at every saved time step.
# for (i, t_val) in enumerate(ts)
#     println("At t = $(t_val):")
#     for point in 1:points
#         println("  Point $(point) spring force: ", norm(spring[i][:, point]))
#     end
# end
# # # Print the drag (drag_force) values for each segment at every saved time step
# # drag_sol = sol[half_drag_force, :]
# # for (i, t_val) in enumerate(ts)
# #     println("At t = $(t_val):")
# #     for seg in 1:segments
# #          # Extract the 3D drag vector for this segment at this time step
# #          println("  Segment $(seg) drag: ", drag_sol[i][:, seg])
# #     end
# # end

# ground_force_spring=(sol[spring_force,10])
# ground_force_drag=(sol[half_drag_force,10])
# ground_force = ground_force_spring + ground_force_drag
# println("Ground force: ", norm.(eachcol(ground_force)))
# Extracting the Apparent Velocity
# -----------------------------
# Note: v_apparent is defined on the segments (each connecting two points).
# To inspect its value over time, we extract it from the solution.
v_app_sol = sol[v_apparent, :]

println("\nApparent velocity (norm) for each segment at each saved time step:")
for (i, t_val) in enumerate(ts)
    println("At t = $(t_val):")
    for seg in 1:segments
         # Compute the norm of the apparent velocity vector for segment 'seg'
         app_norm = norm(v_app_sol[i][:, seg])
         println("  Segment $(seg): norm = $(app_norm)")
    end
end

# Plotting
# Extract solution positions
pos_sol = sol[pos, :]
# Initialize empty lists for all time steps
X, Y, Z = [], [], []
# Extract position data for all points over time
for i in 1:length(pos_sol)
    x, y, z = [], [], []  # temporary lists to store data for current time step
    # Extract P1 (First Point)
    push!(x, pos_sol[i][1])
    push!(y, pos_sol[i][2])
    push!(z, pos_sol[i][3])
    # Extract all other points (P2 to Pn)
    for j in 2:points
       push!(x, pos_sol[i][3*(j-1) + 1])
       push!(y, pos_sol[i][3*(j-1) + 2])
       push!(z, pos_sol[i][3*(j-1) + 3])
    end
    push!(X, x)
    push!(Y, y)
    push!(Z, z)
end

# Animation Setup
lines, sc = nothing, nothing
ylim = (minimum(vcat(Y...))-2, maximum(vcat(Y...))+2)
zlim = (minimum(vcat(Z...))-2, maximum(vcat(Z...))+2)

for i in 1:length(Z)
    z = Z[i]
    y = Y[i]
    global lines, sc
    lines, sc = plot_kite(y, z, ylim, zlim, lines, sc, conn)
    plt.pause(0.01)
    plt.show(block=false)
end
nothing  # Ensures script runs to completion