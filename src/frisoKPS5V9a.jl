#naming and comments made more clear
#different dampings, relative dampings added, Unit damping added

using Timers
tic()
using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra, Timers, Parameters, ControlPlots
using ModelingToolkit: t_nounits as t, D_nounits as D
toc()
include("videoKPS5.jl")

@with_kw mutable struct Settings @deftype Float64
    g_earth::Vector{Float64} = [0.0, 0.0, -9.81]         # gravitational acceleration [m/s²]
    v_wind_tether::Vector{Float64} = [0.0, 15, 0.0]      # wind velocity [m/s]
    v_ro::Float64 = 1.0                                  # reel-out speed [m/s]
    tethersegments::Int64 = 6                            # number of tether segments [-]
    segments::Int64 = 9 + tethersegments                 # total segments [-]
    points::Int64 = 5 + tethersegments                   # total points [-]
    duration::Float64 = 10                               # simulation duration [s]
    save::Bool = false                                   # save animation frames
    K1::Float64 = 614600.0                               # TETHER unit spring constant [N/m]
    K2::Float64 = 10000.0                                # BRIDLE unit spring constant [N/m]
    K3::Float64 = 60000.0                                # KITE unit spring constant [N/m]
    m_kite::Float64 = 6.2 / 4                            # mass of kite [kg]
    m_kcu::Float64 = 8.4                                 # mass of KCU  (at bridle point)[kg]
    rho_tether::Float64 = 724.0                          # density of tether [kg/m³]
    damping::Float64 = 0.9                               # damping coefficient
    rho::Float64 = 1.225                                 # air density [kg/m³]
    cd_tether::Float64 = 0.958                           # drag coefficient of tether
    d_tether::Float64 = 0.004                            # tether diameter [m]
    Cl::Float64 = 0.1                                    # lift coefficient
    S::Float64 = 7.0                                     # surface area [m²]
    kcu_cd::Float64 = 0.47                               # KCU drag coefficient
    kcu_diameter::Float64 = 0.38                         # KCU diameter [m]
end

function calc_initial_state(se)               #first point bridle system, sixth point origin
    POS0 = [0.000   1    -1    0    0    0.000;
            0.000   0     0    2   -2    0.000;
            10.000  13   13   12   12    0.000]
    if se.tethersegments > 1                  # adding extra nodes for tether segments
        extra_nodes = [POS0[:,6] + (POS0[:,1] - POS0[:,6]) * i / se.tethersegments for i in 1:(se.tethersegments-1)]
        POS0 = hcat(POS0, extra_nodes...)
    end     
    VEL0 = zeros(3, se.points)
    POS0, VEL0
end

function model(se)
    POS0, VEL0 = calc_initial_state(se)
    @parameters K1=se.K1 K2=se.K2 K3=se.K3 m_kite=se.m_kite m_kcu=se.m_kcu rho_tether=se.rho_tether 
    @parameters l1=sqrt(10) l2=2.0 l3=sqrt(10) l4=sqrt(8) l5=sqrt(6) l6=sqrt(6) l7=sqrt(8) l8=sqrt(6) l9=sqrt(6) l10=10
    @parameters damping=se.damping
    @parameters rho=se.rho cd_tether=se.cd_tether d_tether=se.d_tether Cl=se.Cl S=se.S 
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
    @variables c_spring(t)[1:se.segments]
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
    eqs2 = vcat(eqs2, acc[:,6] .~ [0.0, 0.0, 0.0])
    # -----------------------------
    # defining the connections and their respective rest lengths, unit spring constants, damping and masses
    # -----------------------------
                                       # connections
    conn = [(1,2), (2,3), (3,1), (1,4), (2,4), (3,4), (1,5), (2,5), (3,5)]
    conn = vcat(conn, [(6+i, 6+i+1) for i in 0:(se.tethersegments-2)]...)
    conn = vcat(conn, [(6+se.tethersegments-1, 1)])
                                       # rest lengths
    rest_lengths = [l1, l2, l3, l4, l5, l6, l7, l8, l9]
    rest_lengths = vcat(rest_lengths, [(l10+se.v_ro*t)/se.tethersegments for _ in 1:se.tethersegments]...)
                                       # unit spring constants (K1 tether, K2 bridle, K3 kite)
    k_segments = [K2, K3, K2, K2, K3, K3, K2, K3, K3]
    k_segments = vcat(k_segments, [K1 for _ in 1:se.tethersegments]...)
                                       # unit damping
                                       # masses
    mass_tether = (d_tether^2)*pi*rho_tether*l10
    mass_tetherpoints = mass_tether/(se.tethersegments+1)
    PointMasses = [se.m_kcu+mass_tetherpoints, se.m_kite, se.m_kite, se.m_kite, se.m_kite]
    PointMasses = vcat(PointMasses, [mass_tetherpoints for _ in 1:se.tethersegments]...)
# -----------------------------
# Equations for Each Segment (Spring Forces, Drag, etc.)
# -----------------------------
    for i in 1:se.segments  
        local eqs = [
           segment[:, i]      ~ pos[:, conn[i][2]] - pos[:, conn[i][1]],
           norm1[i]           ~ norm(segment[:, i]),
           unit_vector[:, i]  ~ -segment[:, i] / norm1[i],
           rel_vel[:, i]      ~ vel[:, conn[i][2]] - vel[:, conn[i][1]],
           spring_vel[i]      ~ -unit_vector[:, i] ⋅ rel_vel[:, i],
           c_spring[i]        ~ (k_segments[i]/rest_lengths[i]) * (0.1 + 0.9*(norm1[i] > rest_lengths[i])),
           spring_force[:, i] ~ (c_spring[i]*(norm1[i] - rest_lengths[i]) + damping * spring_vel[i]) * unit_vector[:, i],
           v_apparent[:, i]   ~ se.v_wind_tether .- (vel[:, conn[i][1]] + vel[:, conn[i][2]]) / 2,
           v_app_perp[:, i]   ~ v_apparent[:, i] - (v_apparent[:, i] ⋅ unit_vector[:, i]) .* unit_vector[:, i],
           norm_v_app[i]      ~ norm(v_app_perp[:, i]),
           half_drag_force[:, i] ~ 0.25 * rho * cd_tether * norm_v_app[i] * (rest_lengths[i]*d_tether) * v_app_perp[:, i]
        ]
        eqs2 = vcat(eqs2, reduce(vcat, eqs))
    end
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
        elseif i != 6                      
            push!(eqs, total_force[:, i] ~ force)
        end
        push!(eqs, acc[:, i] ~ se.g_earth + total_force[:, i] / PointMasses[i])
        eqs2 = vcat(eqs2, reduce(vcat, eqs))
    end

    for i in 1:se.points
        eqs2 = vcat(eqs2, v_app_point[:, i] ~ se.v_wind_tether - vel[:, i])
    end

# Build the ODE System
    @named sys = ODESystem(reduce(vcat, Symbolics.scalarize.(eqs2)), t)
    simple_sys = structural_simplify(sys) 
    simple_sys, pos, vel
end

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

    # Debugging: Print Y and Z arrays
    println("Y coordinates:")
    println(Y)
    println("Z coordinates:")
    println(Z)

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
function main()
    se = Settings()
    simple_sys, pos, vel = model(se)
    sol, elapsed_time = simulate(se, simple_sys, pos, vel)

    # Define connectivity
    conn = [(1,2), (2,3), (3,1), (1,4), (2,4), (3,4), (1,5), (2,5), (3,5)]
    conn = vcat(conn, [(6+i, 6+i+1) for i in 0:(se.tethersegments-2)]...)
    conn = vcat(conn, [(6+se.tethersegments-1, 1)])

    play(se, sol, pos, conn)
    println("Elapsed time: $(elapsed_time) s, speed: $(round(se.duration/elapsed_time)) times real-time")
end

main()
