using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra, Timers, Parameters, ControlPlots
using ModelingToolkit: t_nounits as t, D_nounits as D
using ControlPlots

# -----------------------------
# Simulation Settings
# -----------------------------
@with_kw mutable struct Settings
    g_earth::Vector{Float64} = [0.0, 0.0, -9.81]  # gravitational acceleration [m/s²]
    v_wind_tether::Vector{Float64} = [0.0, 15, 0.0]  # Wind velocity [m/s]
    v_ro = 0.5                                       # reel-out velocity [m/s]
    tethersegments::Int64 = 3                        # number of tether segments [-]
    segments::Int64 = 9 + tethersegments
    points::Int64 = 5 + tethersegments
    duration::Float64 = 10.0                         # simulation time [s]
    save::Bool = false                               # save png files in folder video
    K1::Float64 = 614600                            # spring constant [N]
    K2::Float64 = 10000
    K3::Float64 = 60000
    m_kite::Float64 = 6.2 / 4                       # kite mass [kg]
    m_bridle::Float64 = 8.4                         # bridle mass [kg]
    rho_tether::Float64 = 724                       # density of Dyneema [kg/m³]
    l1::Float64 = sqrt(10)
    l2::Float64 = 2.0
    l3::Float64 = sqrt(10)
    l4::Float64 = sqrt(8)
    l5::Float64 = sqrt(6)
    l6::Float64 = sqrt(6)
    l7::Float64 = sqrt(8)
    l8::Float64 = sqrt(6)
    l9::Float64 = sqrt(6)
    l10::Float64 = 10
    damping::Float64 = 0.9
    rho::Float64 = 1.225
    cd_tether::Float64 = 0.958
    d_tether::Float64 = 0.004
    Cl::Float64 = 0.1
    S::Float64 = 7
end

function calc_initial_state(se)
    # Original positions
    POS0 = [0.000   1    -1    0    0    0.000;
            0.000   0     0    2   -2    0.000;
            10.000  13   13   12   12    0.000]
    if se.tethersegments > 1
        extra_nodes = [POS0[:,6] + (POS0[:,1] - POS0[:,6]) * i / se.tethersegments for i in 1:(se.tethersegments-1)]
        POS0 = hcat(POS0, extra_nodes...)
    end     
    VEL0 = zeros(3, se.points)
    POS0, VEL0
end

function model(se, POS0, VEL0)
    # -----------------------------
    @parameters K1=se.K1  K2=se.K2 K3=se.K3 m_kite=se.m_kite m_bridle=se.m_bridle rho_tether=se.rho_tether
    @parameters l1=se.l1 l2=se.l2 l3=se.l3 l4=se.l4 l5=se.l5 l6=se.l6 l7=se.l7 l8=se.l8 l9=se.l9 l10=se.l10 damping=se.damping 
    @parameters rho=se.rho cd_tether=se.cd_tether d_tether=se.d_tether Cl=se.Cl S=se.S
    # State variables for points
        # State variables for points
    @variables pos(t)[1:3, 1:se.points] = POS0
    @variables vel(t)[1:3, 1:se.points] = VEL0
    @variables acc(t)[1:3, 1:se.points]

    # Define additional variable: apparent velocity for each point.
    @variables v_app_point(t)[1:3, 1:se.points]
    # (v_app_point = v_wind_tether - vel will be added later)

    # Variables for segments (connecting points)
    @variables segment(t)[1:3, 1:se.segments]
    @variables unit_vector(t)[1:3, 1:se.segments]
    @variables norm1(t)[1:se.segments]
    @variables rel_vel(t)[1:3, 1:se.segments]
    @variables spring_vel(t)[1:se.segments]
    @variables c_spring(t)[1:se.segments]
    @variables spring_force(t)[1:3, 1:se.segments]
    @variables v_apparent(t)[1:3, 1:se.segments]
    @variables v_app_perp(t)[1:3, 1:se.segments]
    @variables norm_v_app(t)[1:se.segments]
    @variables half_drag_force(t)[1:3, 1:se.segments]
    @variables total_force(t)[1:3, 1:se.points]

    # -----------------------------
    # Basic Differential Equations for Points
    # -----------------------------
    eqs = vcat(D.(pos) .~ vel,
            D.(vel) .~ acc)
    eqs2 = vcat(eqs...)   # define eqs2 globally
    # Fix point 6 (anchor) acceleration to zero:
    eqs2 = vcat(eqs2, acc[:,6] .~ [0.0, 0.0, 0.0])

    # Connectivity, Rest Lengths, and Stiffness Constants
    conn = [(1,2), (2,3), (3,1), (1,4), (2,4), (3,4), (1,5), (2,5), (3,5)]
    conn = vcat(conn, [(6+i, 6+i+1) for i in 0:(se.tethersegments-2)]...)  # tether segments
    conn = vcat(conn, [(6+se.tethersegments-1, 1)])
    rest_lengths = [l1, l2, l3, l4, l5, l6, l7, l8, l9]
    rest_lengths = vcat(rest_lengths, [(l10 + se.v_ro * t) / se.tethersegments for _ in 1:se.tethersegments]...)
    k_segments = [K2, K3, K2, K2, K3, K3, K2, K3, K3]
    k_segments = vcat(k_segments, [K1 for _ in 1:se.tethersegments]...)

    # Mass Distribution
    mass_tether = (d_tether^2) * π * rho_tether * l10
    mass_tetherpoints = mass_tether / (se.tethersegments + 1)
    PointMasses = [se.m_bridle + mass_tetherpoints, se.m_kite, se.m_kite, se.m_kite, se.m_kite]
    PointMasses = vcat(PointMasses, [mass_tetherpoints for _ in 1:se.tethersegments]...)

    # Define Forces and Update Equations
    for i in 1:se.segments  
        eqs = [
           segment[:, i]      ~ pos[:, conn[i][2]] - pos[:, conn[i][1]],
           norm1[i]           ~ norm(segment[:, i]),
           unit_vector[:, i]  ~ -segment[:, i] / norm1[i],
           rel_vel[:, i]      ~ vel[:, conn[i][2]] - vel[:, conn[i][1]],
           spring_vel[i]      ~ -unit_vector[:, i] ⋅ rel_vel[:, i],
           c_spring[i]        ~ (k_segments[i]/rest_lengths[i]) * (0.01 + 0.99*(norm1[i] > rest_lengths[i])),
           spring_force[:, i] ~ (c_spring[i]*(norm1[i] - rest_lengths[i]) + damping * spring_vel[i]) * unit_vector[:, i],
           v_apparent[:, i]   ~ se.v_wind_tether .- (vel[:, conn[i][1]] + vel[:, conn[i][2]]) / 2,
           v_app_perp[:, i]   ~ v_apparent[:, i] - (v_apparent[:, i] ⋅ unit_vector[:, i]) .* unit_vector[:, i],
           norm_v_app[i]      ~ norm(v_app_perp[:, i]),
           half_drag_force[:, i] ~ 0.25 * rho * cd_tether * norm_v_app[i] * (rest_lengths[i]*d_tether) * v_app_perp[:, i]
        ]
        eqs2 = vcat(eqs2, reduce(vcat, eqs))
    end
    for i in 1:se.points  
        eqs = []  
        # Compute total force acting on each point from spring forces and drag forces.
        force = sum([spring_force[:, j] for j in 1:se.segments if conn[j][2] == i]; init=zeros(3)) - 
                sum([spring_force[:, j] for j in 1:se.segments if conn[j][1] == i]; init=zeros(3)) + 
                sum([half_drag_force[:, j] for j in 1:se.segments if conn[j][1] == i]; init=zeros(3)) + 
                sum([half_drag_force[:, j] for j in 1:se.segments if conn[j][2] == i]; init=zeros(3))
        v_app_point[:, i] ~ se.v_wind_tether - vel[:, i]
        if i in 2:5
            #L = 0.5 * se.rho * se.Cl * se.S * (v_app_point[1, i]^2 + v_app_point[2, i]^2 + v_app_point[3, i]^2)
            L=345
            push!(eqs, total_force[:, i] ~ force + [0.0, 0.0, L]) 
        elseif i != 6        
            push!(eqs, total_force[:, i] ~ force)
        end
        push!(eqs, acc[:, i] ~ se.g_earth + total_force[:, i] / PointMasses[i])
        eqs2 = vcat(eqs2, reduce(vcat, eqs))
        eqs2 = vcat(eqs2, v_app_point[:, i] ~ se.v_wind_tether - vel[:, i])
    end
    
    eqs2 = vcat(eqs2, reduce(vcat, eqs)) 
    @named sys = ODESystem(reduce(vcat, Symbolics.scalarize.(eqs2)), t)
    simple_sys = structural_simplify(sys)
    simple_sys, pos, vel
end

# # Main Simulation Function
# function simulate(se, simple_sys)
#     dt = 0.02
#     tol = 1e-6
#     tspan = (0.0, se.duration)
#     ts    = 0:dt:se.duration
#     prob = ODEProblem(simple_sys, nothing, tspan)
#     elapsed_time = @elapsed sol = solve(prob, Rodas5(); dt=dt, abstol=tol, reltol=tol, saveat=ts)
#     sol, elapsed_time
# end

# # Plotting Function (for animation)
# function play(se, sol)
#     # Add your plotting code here
#     nothing
# end

# Run the simulation
function main()
    se = Settings()
    POS0, VEL0 = calc_initial_state(se)
    eqs = model(se, POS0, VEL0)
    sol = simulate(se, eqs)
    #play(se, sol)
    println("Simulation complete")
end

main()
