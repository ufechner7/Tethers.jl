using Timers
tic()
using LinearAlgebra, ModelingToolkit, OrdinaryDiffEq, ControlPlots 
using ModelingToolkit: t_nounits as t, D_nounits as D
import Plots
toc()
include("video.jl")
# 🔹 Define struct for simulation settings
G_EARTH::Vector{Float64} = [0.0, -9.81] 
F2::Vector{Float64} = [0.0, 13]
F3::Vector{Float64} = [0.0, 13]
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

# basic differential equations
eqs1 = vcat(D.(pos) .~ vel,
            D.(vel) .~ acc,
            )
eqs2 = vcat(eqs1...)
eqs2 = vcat(eqs2, pos[:,1] .~ [0.0,0.0])  
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
    force = sum(
        spring_force[:, j] for j in 1:segments if conn[j][2] == i  # Forces where this point is the end
    ) - sum(
        spring_force[:, j] for j in 1:segments if conn[j][1] == i  # Forces where this point is the start
    )
    # Apply Newton's Second Law: F = ma, ensuring correct force addition
    if i == 2
        push!(eqs, total_force[:, i] ~ force + F2)  # Apply F2 as vector addition
    elseif i == 3
        push!(eqs, total_force[:, i] ~ force + F3)  # Apply F3 as vector addition
    end
    push!(eqs, acc[:, i] ~ G_EARTH + total_force[:, i] / m)  # Compute acceleration with gravity
    # Append the equations for this particle to the global system
    eqs2 = vcat(eqs2, reduce(vcat, eqs))
end

@named sys = ODESystem(reduce(vcat, Symbolics.scalarize.(eqs2)), t)
simple_sys = structural_simplify(sys)
# running the simulation
dt = 0.02
tol = 1e-6
tspan = (0.0, duration)
ts    = 0:dt:duration
prob = ODEProblem(simple_sys, nothing, tspan)
elapsed_time = @elapsed sol = solve(prob, KenCarp4(autodiff=false); dt, abstol=tol, reltol=tol, saveat=ts)
println("Elapsed time: $(elapsed_time) s, speed: $(round(duration/elapsed_time)) times real-time")
#plotting
pos_sol = sol[pos, :]
X=[]
Y=[]
for i in 1:length(pos_sol)
    x2 = pos_sol[i][3]
    x3 = pos_sol[i][5]
    y2 = pos_sol[i][4]
    y3 = pos_sol[i][6]
    x = [0,x2, x3] 
    y = [0,y2, y3] 
    push!(X,x) 
    push!(Y,y)
end
lines,sc = nothing, nothing
xlim=(-5,5)
ylim=(0,10)
for i in 1:length(X)
    x=X[i]
    y=Y[i]
    global lines, sc
    println(x)
    println(y)
    lines,sc = plot_kite(x,y,xlim,ylim,lines,sc)
    #plt.tight_layout()
    plt.pause(0.01)
    plt.show(block=false)
    #sleep(0.05)
 end
nothing