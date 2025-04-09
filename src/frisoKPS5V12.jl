using Timers
tic()
using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra, Timers, Parameters, ControlPlots
using ModelingToolkit: Symbolics, @register_symbolic
using Dierckx
using ModelingToolkit: t_nounits as t, D_nounits as D
toc()
include("videoKPS5.jl")
# -----------------------------
# Interpolating polars using Dierckx
# -----------------------------
alpha_cl= [-180.0, -160.0, -90.0, -20.0, -10.0,  -5.0,  0.0, 20.0, 40.0, 90.0, 160.0, 180.0]
cl_list = [   0.0,    0.5,   0.0,  0.08, 0.125,  0.15,  0.2,  1.0,  1.0,  0.0,  -0.5,   0.0]
alpha_cd = [-180.0, -170.0, -140.0, -90.0, -20.0, 0.0, 20.0, 90.0, 140.0, 170.0, 180.0]
cd_list = [   0.5,    0.5,    0.5,   1.0,   0.2, 0.1,  0.2,  1.0,   0.5,   0.5,   0.5]
function cl_interp(alpha)
    cl_spline = Spline1D(alpha_cl, cl_list)
    return cl_spline(alpha)
end
function cd_interp(alpha)
    cd_spline = Spline1D(alpha_cd, cd_list)
    return cd_spline(alpha)
end

@register_symbolic cl_interp(alpha)
@register_symbolic cd_interp(alpha)
# -----------------------------
@with_kw mutable struct Settings @deftype Float64
    g_earth::Vector{Float64} = [0.0, 0.0, -9.81]         # gravitational acceleration [m/s²]
    v_wind_tether::Vector{Float64} = [13.5, 0.0, 0.0]      # wind velocity [m/s]
    rho::Float64 = 1.225                                 # air density [kg/m³]
    duration::Float64 = 10                           # simulation duration [s]
    save::Bool = false                                   # save animation frames
    # kite
    m_kite::Float64 = 10.58                               # mass of kite [kg]
    S::Float64 = 20.36                                      # surface area [m²]
    kite_width::Float64 =8.16                                      # width of kite [m]
    kite_height::Float64 =3.15                                    # height of kite [m]
    chord_length::Float64 =2.0                                  # chord length [m]
    # KCU + bridle
    bridle_height = 4.9                                  # height of bridle [m]
    d_bridleline::Float64 = 0.00375                                     # bridle line diameter [mm]
    l_bridle::Float64 = 33.4                           # sum of the lengths of the bridle lines [m]
    kcu_cd::Float64 = 0.47                               # KCU drag coefficient
    kcu_diameter::Float64 = 0.38                         # KCU diameter [m]
    m_kcu::Float64 = 11                                 # mass of KCU  (at bridle point)[kg]
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
    v_ro::Float64 = 0.0                                  # reel-out speed [m/s]
    rho_tether::Float64 = 724.0                          # density of tether [kg/m³]
    cd_tether::Float64 = 0.958                           # drag coefficient of tether
    d_tether::Float64 = 0.004                            # tether diameter [m]
    beta::Float64 = 0.296705973      # 17degrees as in uwes sim,  change this name                         # angle XZ plane, between origin and bridle point [rad]
    l_totaltether::Float64 = 40.0                        # tether length [m]
    tethersegments::Int64 = 14                            # number of tether segments [-]
    segments::Int64 = 9 + tethersegments                 # total segments [-]
    points::Int64 = 5 + tethersegments                   # total points [-]
end

# Function to compute apparent velocity in xz-plane and calculate Alpha1p
function compute_alpha1p(v_a, e_z, e_x)
    # Calculate the angle Alpha1p
    v_a_z = v_a ⋅ e_z  # Use the symbolic dot product
    v_a_x = v_a ⋅ e_x  # Use the symbolic dot product
    
    # Use atan for the symbolic computation
    alpha1p = atan(v_a_z, v_a_x)
    return alpha1p
end
function xz_distance(p1, p2)
    return sqrt((p2[1] - p1[1])^2 + (p2[3] - p1[3])^2)
end
function get_kite_points(se)
    # Original kite points in local reference frame
    kitepos0 =                         # KITE points  
    # P1 Bridle        P2                              P3                            P4                            P5
    [0.000         se.chord_length/2             -se.chord_length/2                  0                              0;
    0.000               0                               0                        -se.kite_width/2           se.kite_width/2;
    0.000     se.kite_height+se.bridle_height  se.kite_height+se.bridle_height   se.bridle_height          se.bridle_height]

    beta = se.beta  # Assuming Beta is stored in se.beta
    Y_r = [cos(beta) 0 sin(beta);
                 0    1       0;
          -sin(beta) 0 cos(beta)]    
    # Apply rotation to all points
    kitepos0rot =  Y_r * kitepos0 
    
    return kitepos0rot
end
# Local Kite Reference Frame
function local_kite_reference_frame(P2, P3, P4, P5)
    X = P2 - P3  # X axis is the vector from P3 to P2
    Y = P5 - P4  # Y axis is the vector from P4 to P5
    Z = cross(X, Y)  # Cross product to get Z axis
    e_x = X / norm(X)  # Unit vector along the X axis
    e_y = Y / norm(Y)  # Unit vector along the Y axis
    e_z = Z / norm(Z)  # Unit vector along the Z axis
    return e_x, e_y, e_z
end

# Calculate Initial State
function calc_initial_state(se)  
    p1location = [se.l_totaltether*sin(se.beta) 0 se.l_totaltether*cos(se.beta)]
    kitepos0rot = get_kite_points(se)
    POS0 = kitepos0rot .+ p1location'
    POS0 = hcat(POS0, zeros(3, 1))
    if se.tethersegments > 1
        extra_nodes = [POS0[:,6] + (POS0[:,1] - POS0[:,6]) * i / se.tethersegments for i in 1:(se.tethersegments-1)]
        POS0 = hcat(POS0, extra_nodes...)
    end     
    VEL0 = zeros(3, se.points)
    return POS0, VEL0
end

# Calculate Rest Lengths
function calc_rest_lengths(se)
    POS0, VEL0  = calc_initial_state(se)  
    lengths = [norm(POS0[:,se.conn[i][2]] - POS0[:,se.conn[i][1]]) for i in 1:9]
    l10 = norm(POS0[:,1] - POS0[:,6])
    lengths = vcat(lengths, [(l10+se.v_ro*t)/se.tethersegments for _ in 1:se.tethersegments]...)
    return lengths, l10
end
# Define a function to create reference frame update equations
function create_reference_frame_equations(pos, e_x, e_y, e_z)
    # Calculate vectors for the reference frame
    X = pos[:, 2] - pos[:, 3]  # Vector from P3 to P2
    Y = pos[:, 5] - pos[:, 4]  # Vector from P4 to P5
    Z = cross(X, Y)            # Cross product for Z axis
    # Normalize these vectors to get unit vectors
    norm_X = sqrt(sum(X .^ 2))
    norm_Y = sqrt(sum(Y .^ 2))
    norm_Z = sqrt(sum(Z .^ 2))
    # Create equations to update the reference frame
    ref_frame_eqs = [
        e_x[1] ~ X[1] / norm_X,
        e_x[2] ~ X[2] / norm_X,
        e_x[3] ~ X[3] / norm_X,
        e_y[1] ~ Y[1] / norm_Y,
        e_y[2] ~ Y[2] / norm_Y,
        e_y[3] ~ Y[3] / norm_Y,
        e_z[1] ~ Z[1] / norm_Z,
        e_z[2] ~ Z[2] / norm_Z,
        e_z[3] ~ Z[3] / norm_Z
    ]  
    return ref_frame_eqs
end
# Define the Model
function model(se)
    POS0, VEL0 = calc_initial_state(se)
    rest_lengths, l_tether = calc_rest_lengths(se)
    @parameters K1=se.springconstant_tether K2=se.springconstant_bridle K3=se.springconstant_kite C1=se.damping_tether C2=se.rel_damping_bridle*se.damping_tether C3=se.rel_damping_kite*se.damping_tether
    @parameters m_kite=se.m_kite m_kcu=se.m_kcu rho_tether=se.rho_tether 
    @parameters rho=se.rho cd_tether=se.cd_tether d_tether=se.d_tether S=se.S 
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
    # local kite reference frame
    @variables e_x(t)[1:3]
    @variables e_y(t)[1:3]
    @variables e_z(t)[1:3] 
    @variables alpha1p(t)[1:4]  
    
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
    mass_bridlelines = ((se.d_bridleline/2)^2)*pi*rho_tether*se.l_bridle #total mass entire bridle 
    mass_halfbridleline = mass_bridlelines/8 # half the connection of bridle line to kite (to assign to each kitepoint) so the other 4 halves get assigned to bridlepoint 
    mass_tether = ((d_tether/2)^2)*pi*rho_tether*l_tether
    mass_tetherpoints = mass_tether/(se.tethersegments+1)
    mass_bridlepoint = 4*mass_halfbridleline + se.m_kcu + mass_tetherpoints # 4 bridle connections, kcu and tether
    m_kitepoints = (se.m_kite/4) + mass_halfbridleline 
    PointMasses = [mass_bridlepoint, m_kitepoints, m_kitepoints, m_kitepoints, m_kitepoints]
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
            norm_v_app[i]      ~ norm(v_app_perp[:, i])
        ]
        if i > 9 # tether segments
            push!(eqs, half_drag_force[:, i] ~ 0.25 * rho * cd_tether * norm_v_app[i] * (rest_lengths[i]*d_tether) * v_app_perp[:, i])
        elseif i in [1, 3, 4, 7] # bridle lines, try to find Cd_bridlelines later
            push!(eqs, half_drag_force[:, i] ~ 0.25 * rho * cd_tether * norm_v_app[i] * (rest_lengths[i]*se.d_bridleline) * v_app_perp[:, i])
        else # kite
            push!(eqs, half_drag_force[:, i] ~ zeros(3))
        end
        eqs2 = vcat(eqs2, reduce(vcat, eqs))
    end
    # -----------------------------
    # Reference Frame
    # -----------------------------
    ref_frame_eqs = create_reference_frame_equations(pos, e_x, e_y, e_z)
    eqs2 = vcat(eqs2, ref_frame_eqs)

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
        elseif i in 2:5 #the kite points that get Aero Forces
            v_a = se.v_wind_tether - vel[:, i]           # Apparent wind velocity
            v_app_mag_squared = v_app_point[1, i]^2 + v_app_point[2, i]^2 + v_app_point[3, i]^2
            alpha1p_i = compute_alpha1p(v_a, e_z, e_x)   # Calculate Alpha1p at this time step
            eqs2 = vcat(eqs2, alpha1p[i-1] ~ alpha1p_i)  # Add the equation for Alpha1p for each of 4 kite points (first bering bridle so i-1)   
            # getting Cl and Cd
            Cl = cl_interp(alpha1p_i)            
            Cd = cd_interp(alpha1p_i)

            # Lift calculation
            L_perpoint = (1/4) * 0.5 * rho * Cl * S * (v_app_mag_squared)
            # Cross product and normalization
            cross_vapp_X_e_y = cross(v_app_point[:, i], e_y)
            normcross_vapp_X_e_y = norm(cross_vapp_X_e_y)
            L_direction = cross_vapp_X_e_y / normcross_vapp_X_e_y
            # Final lift force vector
            L = L_perpoint * L_direction

            # Drag calculation
            D_perpoint = (1/4) * 0.5 * rho * Cd * S * v_app_mag_squared
            # Create drag direction components
            D_direction = [v_app_point[1, i] / norm(v_app_point[:, i]), v_app_point[2, i] / norm(v_app_point[:, i]), v_app_point[3, i] / norm(v_app_point[:, i])]
            # Final drag force vector components
            D = D_perpoint * D_direction
            
            # Total aerodynamic force
            Fa = [L[1]+ D[1], L[2]+ D[2], L[3]+ D[3]]
            
            push!(eqs, total_force[1, i] ~ force[1] + Fa[1])
            push!(eqs, total_force[2, i] ~ force[2] + Fa[2])
            push!(eqs, total_force[3, i] ~ force[3] + Fa[3])
        elseif i != 6                      
            push!(eqs, total_force[:, i] ~ force)
        end
        push!(eqs, acc[:, i] ~ se.g_earth + total_force[:, i] / PointMasses[i])
        eqs2 = vcat(eqs2, reduce(vcat, eqs))
        eqs2 = vcat(eqs2, v_app_point[:, i] ~ se.v_wind_tether - vel[:, i])
    end
    @named sys = ODESystem(reduce(vcat, Symbolics.scalarize.(eqs2)), t)
    simple_sys = structural_simplify(sys) 
    simple_sys, pos, vel, conn, e_x, e_y, e_z, v_app_point, alpha1p 
end
# -----------------------------
# Simulation Function
# -----------------------------
function simulate(se, simple_sys, pos, vel, alpha1p; prn=false)
    dt = 0.02
    tol = 1e-6
    tspan = (0.0, se.duration)
    ts = 0:dt:se.duration
    prob = ODEProblem(simple_sys, nothing, tspan)
    elapsed_time = @elapsed sol = solve(prob, Rodas5(autodiff=false); dt=dt, abstol=tol, reltol=tol, saveat=ts)

    # Debugging: Print the solution
    if prn
        println("Solution (sol):")
        #println(sol)
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
    X, Z = [], [] 
    for i in 1:length(pos_sol)
        x, z = [], []  # Renamed from x, y to y, z
        for j in 1:se.points
            push!(x, pos_sol[i][1, j])  # Extract y-coordinate for point j at time step i
            push!(z, pos_sol[i][3, j])  # Extract z-coordinate for point j at time step i
        end
        push!(X, x)
        push!(Z, z)
    end
    lines, sc = nothing, nothing
    zlim = (minimum(vcat(Z...)) - 2, maximum(vcat(Z...)) + 2)  # Renamed from ylim to zlim
    xlim = (minimum(vcat(X...)) - 2, maximum(vcat(X...)) + 2)  # Renamed from xlim to ylim

    for i in 1:length(X)
        x = X[i]
        z = Z[i]
        lines, sc = plot_kite(x, z, xlim, zlim, lines, sc, conn)  # Changed order of arguments if needed
        plt.pause(0.01)
        plt.show(block=false)
    end
    nothing
end
# -----------------------------
# plotting AOA over time
# -----------------------------
function plot_alpha1p_control(sol, alpha1p)
    # Extract time values
    time = sol.t
    
    # Extract alpha1p values for all 4 points
    alpha1p_values = [rad2deg.(sol[alpha1p[i], :]) for i in 1:4]
    
    # Create a figure
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)
    
    # Plot each alpha1p value
    ax.plot(time, alpha1p_values[1], label="Kite Point 1")
    ax.plot(time, alpha1p_values[2], label="Kite Point 2")
    ax.plot(time, alpha1p_values[3], label="Kite Point 3")
    ax.plot(time, alpha1p_values[4], label="Kite Point 4")
    
    # Add labels and title
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Angle of Attack (rad)")
    ax.set_title("Angle of Attack over Time")
    ax.legend()
    ax.grid(true)
    
    # Display the plot
    plt.tight_layout()
    plt.show()
    
    # Save the figure if needed
    plt.savefig("alpha1p_over_time.png")
    
    return fig
end
# Main function to run the simulation
function main()
    se = Settings()
    simple_sys, pos, vel, conn, e_x, e_y, e_z, v_app_point, alpha1p = model(se)
    sol, elapsed_time = simulate(se, simple_sys, pos, vel, alpha1p; prn=true)
    
    println("Alpha1p over time:")
    for i in 1:length(sol.t)
        println("t = $(sol.t[i]): Alpha1p = $(sol[alpha1p, i])")
    end
    plot_alpha1p_control(sol, alpha1p)
    play(se, sol, pos, conn)
    println("Elapsed time: $(elapsed_time) s, speed: $(round(se.duration/elapsed_time)) times real-time")
end

main()