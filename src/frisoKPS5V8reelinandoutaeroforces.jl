#added aero forces and KCU Drag
using Timers
tic()
using LinearAlgebra, ModelingToolkit, OrdinaryDiffEq, ControlPlots 
using ModelingToolkit: t_nounits as t, D_nounits as D
toc()
include("videoKPS5.jl")

# -----------------------------
# Simulation Settings
# -----------------------------
# Coordinate system: [x,y,z] where x is heading, z is up, and y is perpendicular.
G_EARTH = [0.0, 0.0, -9.81]
v_wind_tether = [0.0, 15, 0.0]
V_RO = -0.5
tethersegments = 6
segments = 9 + tethersegments    
points = 5 + tethersegments  
duration = 10              # Simulation time [s]
save = false                 # Whether to save animation frames

# Original positions
POS0 = [0.000   1    -1    0    0    0.000;
        0.000   0     0    2   -2    0.000;
        10.000  13   13   12   12    0.000]
if tethersegments > 1
    extra_nodes = [POS0[:,6] + (POS0[:,1] - POS0[:,6]) * i / tethersegments for i in 1:(tethersegments-1)]
    POS0 = hcat(POS0, extra_nodes...)
end     
VEL0 = zeros(3, points)

# -----------------------------
# Define Model Parameters and Variables
# -----------------------------
@parameters K1=614600  K2=10000 K3=60000 m_kite=6.2/4 m_bridle=8.4 rho_tether=724 
@parameters l1=sqrt(10) l2=2.0 l3=sqrt(10) l4=sqrt(8) l5=sqrt(6) l6=sqrt(6) l7=sqrt(8) l8=sqrt(6) l9=sqrt(6) l10=10
@parameters damping=0.9
@parameters rho=1.225 cd_tether=0.958 d_tether=0.004 Cl=0.1 S=7 
@parameters kcu_cd = 0.47 kcu_diameter = 0.38


# State variables for points
@variables pos(t)[1:3, 1:points] = POS0
@variables vel(t)[1:3, 1:points] = VEL0
@variables acc(t)[1:3, 1:points]

# Define additional variable: apparent velocity for each point.
@variables v_app_point(t)[1:3, 1:points]
# (v_app_point = v_wind_tether - vel will be added later)

# Variables for segments (connecting points)
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
@variables drag_force(t)[1:3, 1:segments]  # (if needed)
@variables total_force(t)[1:3, 1:points]

# -----------------------------
# Basic Differential Equations for Points
# -----------------------------
eqs1 = vcat(D.(pos) .~ vel,
            D.(vel) .~ acc)
global eqs2 = vcat(eqs1...)   # define eqs2 globally
# Fix point 6 (anchor) acceleration to zero:
eqs2 = vcat(eqs2, acc[:,6] .~ [0.0, 0.0, 0.0])

# -----------------------------
# Connectivity, Rest Lengths, and Stiffness Constants
# -----------------------------
conn = [(1,2), (2,3), (3,1), (1,4), (2,4), (3,4), (1,5), (2,5), (3,5)]
conn = vcat(conn, [(6+i, 6+i+1) for i in 0:(tethersegments-2)]...)  # tether segments
conn = vcat(conn, [(6+tethersegments-1, 1)])
rest_lengths = [l1, l2, l3, l4, l5, l6, l7, l8, l9]
rest_lengths = vcat(rest_lengths, [(l10+V_RO*t)/tethersegments for _ in 1:tethersegments]...)
k_segments = [K2, K3, K2, K2, K3, K3, K2, K3, K3]
k_segments = vcat(k_segments, [K1 for _ in 1:tethersegments]...)

# -----------------------------
# Equations for Each Segment (Spring Forces, Drag, etc.)
# -----------------------------
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
       half_drag_force[:, i] ~ 0.25 * rho * cd_tether * norm_v_app[i] * (rest_lengths[i]*d_tether) * v_app_perp[:, i]
    ]
    eqs2 = vcat(eqs2, reduce(vcat, eqs))
end

# -----------------------------
# Mass Distribution
# -----------------------------
mass_tether = (d_tether^2)*pi*rho_tether*l10
mass_tetherpoints = mass_tether/(tethersegments+1)
PointMasses = [m_bridle+mass_tetherpoints, m_kite, m_kite, m_kite, m_kite]
PointMasses = vcat(PointMasses, [mass_tetherpoints for _ in 1:tethersegments]...)

# -----------------------------
# Force Balance at Each Point
# -----------------------------
for i in 1:points  
    global eqs2
    local eqs = []  
    # Compute total force acting on each point from spring forces and drag forces.
    force = sum([spring_force[:, j] for j in 1:segments if conn[j][2] == i]; init=zeros(3)) -
            sum([spring_force[:, j] for j in 1:segments if conn[j][1] == i]; init=zeros(3)) +
            sum([half_drag_force[:, j] for j in 1:segments if conn[j][1] == i]; init=zeros(3)) +
            sum([half_drag_force[:, j] for j in 1:segments if conn[j][2] == i]; init=zeros(3))
    #ExternalForces = [F1, F2, F3, F4, F5]
    v_app_point[:, i] ~ v_wind_tether - vel[:, i]
    if i == 1

        area_kcu = pi * ((kcu_diameter / 2) ^ 2)
        Dx = 0.5*rho*kcu_cd *area_kcu*(v_app_point[1, i]*v_app_point[1, i])
        Dy = 0.5*rho*kcu_cd *area_kcu*(v_app_point[2, i]*v_app_point[2, i])
        Dz = 0.5*rho*kcu_cd *area_kcu*(v_app_point[3, i]*v_app_point[3, i])
        D = [Dx, Dy, Dz]

        push!(eqs, total_force[:, i] ~ force + D)
    elseif i in 2:5
        L=0.5*rho*Cl*S*(v_app_point[1, i]*v_app_point[1, i] + v_app_point[2, i]*v_app_point[2, i] + v_app_point[3, i]*v_app_point[3, i])
        push!(eqs, total_force[:, i] ~ force + [0.0, 0.0, L]) 
    elseif i != 6 && i != 1        # (optional additional constraint)                  
        push!(eqs, total_force[:, i] ~ force)
    end
    push!(eqs, acc[:, i] ~ G_EARTH + total_force[:, i] / PointMasses[i])
    eqs2 = vcat(eqs2, reduce(vcat, eqs))
end

# -----------------------------
# Define Apparent Velocity for Each Point
# -----------------------------
for i in 1:points
    global eqs2
    eqs2 = vcat(eqs2, v_app_point[:, i] ~ v_wind_tether - vel[:, i])
end

# -----------------------------
# Build and Solve the ODE System
# -----------------------------
@named sys = ODESystem(reduce(vcat, Symbolics.scalarize.(eqs2)), t)
simple_sys = structural_simplify(sys) 

dt = 0.02
tol = 1e-6
tspan = (0.0, duration)
ts = 0:dt:duration
prob = ODEProblem(simple_sys, nothing, tspan)
elapsed_time = @elapsed sol = solve(prob, Rodas5(); dt=dt, abstol=tol, reltol=tol, saveat=ts)
println("Elapsed time: $(elapsed_time) s, speed: $(round(duration/elapsed_time)) times real-time")

# -----------------------------
# Extract and Print the Apparent Velocity for Each Point
# -----------------------------
v_app_point_sol = sol[v_app_point, :]

println("\nApparent velocity (norm) for each point at each saved time step:")
for (i, t_val) in enumerate(ts)
    println("At t = $(t_val):")
    for pt in 1:points
         app_norm = norm(v_app_point_sol[i][:, pt])
         println("  Point $(pt): norm = $(app_norm)")
    end
    for pt in 1
        v_app = v_app_point_sol[i][:, pt]
        Dx = 0.032648334 * (v_app[1]^2) 
        Dy = 0.032648334 * (v_app[2]^2)
        Dz = 0.032648334 * (v_app[3]^2)
        D = [Dx, Dy, Dz]
        println("  Point $(pt): D = $(D)")
    end
end

# -----------------------------
# Plotting (Animation)
# -----------------------------
pos_sol = sol[pos, :]
X, Y, Z = [], [], []
for i in 1:length(pos_sol)
    x, y, z = [], [], []
    push!(x, pos_sol[i][1])
    push!(y, pos_sol[i][2])
    push!(z, pos_sol[i][3])
    for j in 2:points
       push!(x, pos_sol[i][3*(j-1) + 1])
       push!(y, pos_sol[i][3*(j-1) + 2])
       push!(z, pos_sol[i][3*(j-1) + 3])
    end
    push!(X, x)
    push!(Y, y)
    push!(Z, z)
end

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
nothing  # Ensures the script runs to completion

