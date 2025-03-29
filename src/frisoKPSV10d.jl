# define 'kite reference frame' and write down in typst.app_norm
# figure out how AOA defined
# assign L and D to it 
# add different plotting ControlPlots
using Timers
tic()
using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra, Timers, Parameters, ControlPlots
using ModelingToolkit: t_nounits as t, D_nounits as D
using Dierckx
toc()
include("videoKPS5.jl")

@with_kw mutable struct Settings @deftype Float64
    g_earth::Vector{Float64} = [0.0, 0.0, -9.81]         # gravitational acceleration [m/s²]
    v_wind_tether::Vector{Float64} = [15, 15, 0.0]      # wind velocity [m/s]
    rho::Float64 = 1.225                                 # air density [kg/m³]
    duration::Float64 = 10                               # simulation duration [s]
    save::Bool = false                                   # save animation frames
    alpha::Float64 = 2.5                                 # angle of attack [rad]
    #kite
    m_kite::Float64 = 6.2                                # mass of kite [kg]
    S::Float64 = 7.0                                     # surface area [m²]
    kite_width = 4                                       # width of kite [m]
    kite_height = 1.0                                    # height of kite [m]
    bridle_height = 2.0                                  # height of bridle [m]
    chord_length = 2.0                                  # chord length [m]
    conn::Vector{Tuple{Int, Int}} = [(1,2), (2,3),
    (3,1), (1,4), (2,4), (3,4), (1,5), (2,5), (3,5)]    # connections between KITE points  
    # spring constants
    springconstant_tether::Float64 = 614600.0            # TETHER unit spring constant [N]
    springconstant_bridle::Float64 = 10000.0             # BRIDLE unit spring constant [N]
    springconstant_kite::Float64 = 60000.0               # KITE unit spring constant [N]
    # damping
    damping_tether::Float64 = 473                        # TETHER unit damping coefficient [Ns]
    rel_damping_kite::Float64 = 6.0                      # KITE relative unit damping coefficient [-]
    rel_damping_bridle:: Float64 = 6.0                   # BRIDLE relative unit damping coefficient [-]
    # tether
    rho_tether::Float64 = 724.0                          # density of tether [kg/m³]
    cd_tether::Float64 = 0.958                           # drag coefficient of tether
    d_tether::Float64 = 0.004                            # tether diameter [m]
    beta::Float64 = pi/8                                 # angle XZ plane, between origin and bridle point [rad]
    l_totaltether::Float64 = 10.0                        # tether length [m]
    v_ro::Float64 = 1.0                                  # reel-out speed [m/s]
    tethersegments::Int64 = 6                            # number of tether segments [-]
    segments::Int64 = 9 + tethersegments                 # total segments [-]
    points::Int64 = 5 + tethersegments                   # total points [-]
    #KCU
    kcu_cd::Float64 = 0.47                               # KCU drag coefficient
    kcu_diameter::Float64 = 0.38                         # KCU diameter [m]
    m_kcu::Float64 = 8.4                                 # mass of KCU  (at bridle point)[kg]
    #polars
    alpha_cl::Vector{Float64} = [-180.0, -160.0, -90.0, -20.0, -10.0,  -5.0,  0.0, 20.0, 40.0, 90.0, 160.0, 180.0]
    cl_list::Vector{Float64}  = [   0.0,    0.5,   0.0,  0.08, 0.125,  0.15,  0.2,  1.0,  1.0,  0.0,  -0.5,   0.0]
    alpha_cd::Vector{Float64} = [-180.0, -170.0, -140.0, -90.0, -20.0, 0.0, 20.0, 90.0, 140.0, 170.0, 180.0]
    cd_list::Vector{Float64}  = [   0.5,    0.5,    0.5,   1.0,   0.2, 0.1,  0.2,  1.0,   0.5,   0.5,   0.5]
end
# ----------------------------
# get kite points
 #-----------------------------
function get_kite_points(se)
    KITEPOS0 =                         # KITE points  
    # P1 Bridle        P2                              P3                            P4                     P5
    [0.000         se.chord_length/2             -se.chord_length/2                  0                              0;
    0.000               0                               0                        se.kite_width/2           -se.kite_width/2;
    0.000     se.kite_height+se.bridle_height  se.kite_height+se.bridle_height   se.bridle_height          se.bridle_height]  
    KITEPOS0
end
# -----------------------------
# Calculate Initial State
# -----------------------------
function calc_initial_state(se)  ##MAKE VARIABLES LOWERCASe
    p1location = [se.l_totaltether*sin(se.beta) 0 se.l_totaltether*cos(se.beta)]
    KITEPOS0 = get_kite_points(se)
    POS0 = KITEPOS0 .+ p1location'
    POS0 = hcat(POS0, zeros(3, 1))
    if se.tethersegments > 1                  # adding extra nodes for tether segments
        extra_nodes = [POS0[:,6] + (POS0[:,1] - POS0[:,6]) * i / se.tethersegments for i in 1:(se.tethersegments-1)]
        POS0 = hcat(POS0, extra_nodes...)
    end     
    VEL0 = zeros(3, se.points)
    POS0, VEL0
end
# -----------------------------
# Calculate Rest Lengths
# -----------------------------
function calc_rest_lengths(se)
    POS0, VEL0 = calc_initial_state(se)  
    # Calculate lengths for the first 9 segments using conn
    lengths = [norm(POS0[:,se.conn[i][2]] - POS0[:,se.conn[i][1]]) for i in 1:9]
    # Calculate l10 separately, the tether lenghts l_tether
    l10 = norm(POS0[:,1] - POS0[:,6])
    # Add tether segment lengths
    lengths = vcat(lengths, [(l10+se.v_ro*t)/se.tethersegments for _ in 1:se.tethersegments]...)
    return lengths, l10
end
# -----------------------------
# Define the Model
# -----------------------------
function model(se)
    POS0, VEL0 = calc_initial_state(se)
    rest_lengths, l_tether = calc_rest_lengths(se)
    # K unit spring constant (K1 tether, K2 bridle, K3 kite); C unit damping constant (C1 tether, C2 bridle, C3 kite) 
    @parameters K1=se.springconstant_tether K2=se.springconstant_bridle K3=se.springconstant_kite C1=se.damping_tether C2=se.rel_damping_bridle*se.damping_tether C3=se.rel_damping_kite*se.damping_tether
    @parameters m_kite=se.m_kite m_kcu=se.m_kcu rho_tether=se.rho_tether 
    @parameters rho=se.rho cd_tether=se.cd_tether d_tether=se.d_tether S=se.S #Cl=se.Cl
    @parameters kcu_cd=se.kcu_cd kcu_diameter=se.kcu_diameter
    @variables pos(t)[1:3, 1:se.points] = POS0
    @variables vel(t)[1:3, 1:se.points] = VEL0
    @variables acc(t)[1:3, 1:se.points]
    @variables v_app_point(t)[1:3, 1:se.points]
    @variables segment(t)[1:3, 1:se.segments]
    @variables unit_vector(t)[1:3, 1:se.segments]
    @variables norm1(t)[1:se.segments]
    @variables rel_vel(t)[1:3, 1:se.segments]
    @variables spring_vel(t)[1:se.segments]
    @variables k_spring(t)[1:se.segments]
    @variables spring_force(t)[1:3, 1:se.segments]
    @variables v_apparent(t)[1:3, 1:se.segments]
    @variables v_app_perp(t)[1:3, 1:se.segments]
    @variables norm_v_app(t)[1:se.segments]
    @variables half_drag_force(t)[1:3, 1:se.segments]
    @variables drag_force(t)[1:3, 1:se.segments]
    @variables total_force(t)[1:3, 1:se.points]

    eqs1 = vcat(D.(pos) .~ vel,
                D.(vel) .~ acc)
    eqs2 = vcat(eqs1...)
    eqs2 = vcat(eqs2, acc[:,6] .~ [0.0, 0.0, 0.0])      # origin is six, make 6 not being hardcoded
    # -----------------------------
    # defining the connections and their respective rest lengths, unit spring constants, damping and masses
    # -----------------------------
                                       # connections   adding segment connections, from origin to bridle 
    conn = vcat(se.conn, [(6+i, 6+i+1) for i in 0:(se.tethersegments-2)]...)
    conn = vcat(conn, [(6+se.tethersegments-1, 1)]) # final connection last tether point to bridle point
                                       # unit spring constants (K1 tether, K2 bridle, K3 kite)
    k_segments = [K2, K3, K2, K2, K3, K3, K2, K3, K3]
    k_segments = vcat(k_segments, [K1 for _ in 1:se.tethersegments]...)
                                       # unit damping constants (C1 tether, C2 bridle, C3 kite)
    c_segments = [C2, C3, C2, C2, C3, C3, C2, C3, C3]
    c_segments = vcat(c_segments, [C1 for _ in 1:se.tethersegments]...)
                                       # masses
    mass_tether = (d_tether^2)*pi*rho_tether*l_tether
    mass_tetherpoints = mass_tether/(se.tethersegments+1)
    m_kitepoints = se.m_kite/4
    PointMasses = [se.m_kcu+mass_tetherpoints, m_kitepoints, m_kitepoints, m_kitepoints, m_kitepoints]
    PointMasses = vcat(PointMasses, [mass_tetherpoints for _ in 1:se.tethersegments]...)
# -----------------------------
# Equations for Each Segment (Spring Forces, Drag, etc.)
# -----------------------------
    for i in 1:se.segments  
        eqs = [
           segment[:, i]      ~ pos[:, conn[i][2]] - pos[:, conn[i][1]],
           norm1[i]           ~ norm(segment[:, i]),
           unit_vector[:, i]  ~ -segment[:, i] / norm1[i],
           rel_vel[:, i]      ~ vel[:, conn[i][2]] - vel[:, conn[i][1]],
           spring_vel[i]      ~ -unit_vector[:, i] ⋅ rel_vel[:, i],
           k_spring[i]        ~ (k_segments[i]/rest_lengths[i]) * (0.1 + 0.9*(norm1[i] > rest_lengths[i])),
           spring_force[:, i] ~ (k_spring[i]*(norm1[i] - rest_lengths[i]) + c_segments[i] * spring_vel[i]) * unit_vector[:, i],
           v_apparent[:, i]   ~ se.v_wind_tether .- (vel[:, conn[i][1]] + vel[:, conn[i][2]]) / 2,
           v_app_perp[:, i]   ~ v_apparent[:, i] - (v_apparent[:, i] ⋅ unit_vector[:, i]) .* unit_vector[:, i],
           norm_v_app[i]      ~ norm(v_app_perp[:, i]),
           half_drag_force[:, i] ~ 0.25 * rho * cd_tether * norm_v_app[i] * (rest_lengths[i]*d_tether) * v_app_perp[:, i]
        ]
        eqs2 = vcat(eqs2, reduce(vcat, eqs))
    end
    # -----------------------------
    # getting Cl and Cd
    # -----------------------------
    alpha_cl=se.alpha_cl
    cl_list = se.cl_list  
    cl_spline = Spline1D(alpha_cl, cl_list)
    Cl = cl_spline(se.alpha)
    alpha_cd = se.alpha_cd
    cd_list = se.cd_list
    cd_spline = Spline1D(alpha_cd, cd_list)
    Cd = cd_spline(se.alpha)
    # -----------------------------
    # Force Balance at Each Point
    # -----------------------------
    for i in 1:se.points  
        eqs = []  
        force = sum([spring_force[:, j] for j in 1:se.segments if conn[j][2] == i]; init=zeros(3)) -
                sum([spring_force[:, j] for j in 1:se.segments if conn[j][1] == i]; init=zeros(3)) +
                sum([half_drag_force[:, j] for j in 1:se.segments if conn[j][1] == i]; init=zeros(3)) +
                sum([half_drag_force[:, j] for j in 1:se.segments if conn[j][2] == i]; init=zeros(3))
        v_app_point[:, i] ~ se.v_wind_tether - vel[:, i]
        if i == 1                                    #KCU drag at bridle point
            area_kcu = pi * ((kcu_diameter / 2) ^ 2)
            Dx_kcu = 0.5*rho*kcu_cd *area_kcu*(v_app_point[1, i]*v_app_point[1, i])
            Dy_kcu = 0.5*rho*kcu_cd *area_kcu*(v_app_point[2, i]*v_app_point[2, i])
            Dz_kcu = 0.5*rho*kcu_cd *area_kcu*(v_app_point[3, i]*v_app_point[3, i])
            D = [Dx_kcu, Dy_kcu, Dz_kcu]
            push!(eqs, total_force[:, i] ~ force + D)
        elseif i in 2:5
            L_perpoint=(1/4)*0.5*rho*Cl*S*(v_app_point[1, i]*v_app_point[1, i] + v_app_point[2, i]*v_app_point[2, i] + v_app_point[3, i]*v_app_point[3, i])
            Dx_kite = (1/4)*0.5*rho*Cd*S*(v_app_point[1, i]*v_app_point[1, i])
            Dy_kite = (1/4)*0.5*rho*Cd*S*(v_app_point[2, i]*v_app_point[2, i])
            push!(eqs, total_force[:, i] ~ force + [Dx_kite, Dy_kite, L_perpoint])    
        elseif i != 6                      
            push!(eqs, total_force[:, i] ~ force)
        end
        push!(eqs, acc[:, i] ~ se.g_earth + total_force[:, i] / PointMasses[i])
        eqs2 = vcat(eqs2, reduce(vcat, eqs))
        eqs2 = vcat(eqs2, v_app_point[:, i] ~ se.v_wind_tether - vel[:, i])
    end

# Build the ODE System
    @named sys = ODESystem(reduce(vcat, Symbolics.scalarize.(eqs2)), t)
    simple_sys = structural_simplify(sys) 
    simple_sys, pos, vel, conn
end
# -----------------------------
# Simulation Function
# -----------------------------
function simulate(se, simple_sys, pos, vel; prn=false)
    dt = 0.02
    tol = 1e-6
    tspan = (0.0, se.duration)
    ts = 0:dt:se.duration
    prob = ODEProblem(simple_sys, nothing, tspan)
    elapsed_time = @elapsed sol = solve(prob, Rodas5(); dt=dt, abstol=tol, reltol=tol, saveat=ts)

    # Debugging: Print the solution
    if prn
        println("Solution (sol):")
        println(sol)
        # Debugging: Print positions and velocities at specific time steps
        for (i, t_val) in enumerate(ts)
            if i % 10 == 1  # Print every 10th time step
                println("\nTime = $t_val")
                println("Positions (pos):")
                println(sol[pos, i])
                println("Velocities (vel):")
                println(sol[vel, i])
            end
        end
    end
    sol, elapsed_time
end
# -----------------------------
# Plotting Function (for animation)
# -----------------------------
function play(se, sol, pos, conn)
    pos_sol = sol[pos, :]
    Y, Z = [], []  # Renamed from X, Y to Y, Z
    for i in 1:length(pos_sol)
        y, z = [], []  # Renamed from x, y to y, z
        for j in 1:se.points
            push!(y, pos_sol[i][2, j])  # Extract y-coordinate for point j at time step i
            push!(z, pos_sol[i][3, j])  # Extract z-coordinate for point j at time step i
        end
        push!(Y, y)
        push!(Z, z)
    end
    lines, sc = nothing, nothing
    zlim = (minimum(vcat(Z...)) - 2, maximum(vcat(Z...)) + 2)  # Renamed from ylim to zlim
    ylim = (minimum(vcat(Y...)) - 2, maximum(vcat(Y...)) + 2)  # Renamed from xlim to ylim

    for i in 1:length(Y)
        y = Y[i]
        z = Z[i]
        lines, sc = plot_kite(y, z, ylim, zlim, lines, sc, conn)  # Changed order of arguments if needed
        plt.pause(0.01)
        plt.show(block=false)
    end
    nothing
end
# -----------------------------
# Main function to run the simulation
# -----------------------------
function main()
    se = Settings()
    simple_sys, pos, vel, conn = model(se)
    sol, elapsed_time = simulate(se, simple_sys, pos, vel)
    play(se, sol, pos, conn)
    println("Elapsed time: $(elapsed_time) s, speed: $(round(se.duration/elapsed_time)) times real-time")
end

main()