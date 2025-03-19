using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra, Plots, Parameters

@with_kw mutable struct Settings3Point @deftype Float64
    g::Float64 = 9.81   # Gravity [m/s²]
    m::Float64 = 1.0    # Mass of P2 and P3 [kg]
    k::Float64 = 1.0    # Spring constant [N/m]
    L0::Float64 = 2.0   # Rest length of springs [m]
    F2::Float64 = 10.0  # Vertical force on P2 [N]
    F3::Float64 = 10.0  # Vertical force on P3 [N]
    duration::Float64 = 10.0  # Simulation duration [s]
end

function initial_conditions(se)
    P1 = [0.0, 0.0]
    P2 = [-1.0, 2.0]
    P3 = [1.0, 2.0]
    V0 = zeros(6)  # Initial velocities [0, 0, 0, 0, 0, 0]
    return vcat(P2, P3, V0)  # [x2, y2, x3, y3, vx2, vy2, vx3, vy3]
end

function equations_of_motion!(du, u, se, t)
    x2, y2, x3, y3, vx2, vy2, vx3, vy3 = u
    P1 = [0.0, 0.0]  # Fixed mass point
    P2 = [x2, y2]
    P3 = [x3, y3]
    
    function spring_force(p1, p2)
        r = p2 - p1
        d = norm(r)
        return -se.k * (d - se.L0) * (r / d)  # Hooke's law
    end
    
    F_s1 = spring_force(P1, P2)  # Force from P1 to P2
    F_s2 = spring_force(P2, P3)  # Force from P2 to P3
    F_s3 = spring_force(P3, P1)  # Force from P3 to P1
    
    # Net force on P2 and P3
    F_net2 = F_s2 - F_s1 + [0, se.F2]  # Adding external vertical force F2
    F_net3 = F_s3 - F_s2 + [0, se.F3]  # Adding external vertical force F3
    
    # Newton's second law: F = m * a => a = F / m
    a2 = F_net2 / se.m
    a3 = F_net3 / se.m
    
    du[1:2] .= vx2, vy2  # dx2/dt, dy2/dt
    du[3:4] .= vx3, vy3  # dx3/dt, dy3/dt
    du[5:6] .= a2  # dvx2/dt, dvy2/dt
    du[7:8] .= a3  # dvx3/dt, dvy3/dt
end

function simulate_system(se)
    u0 = initial_conditions(se)
    tspan = (0.0, se.duration)
    prob = ODEProblem(equations_of_motion!, u0, tspan, se)
    sol = solve(prob, Tsit5(), saveat=0.02)
    return sol
end

function animate_system(se, sol)
    anim = @animate for i in eachindex(sol.t)
        x2, y2, x3, y3 = sol[1:4, i]
        scatter!([0, x2, x3], [0, y2, y3], markersize=5, label="Masses", color=:blue)
        plot!([0, x2], [0, y2], label="Spring S1", linewidth=2, color=:red)
        plot!([x2, x3], [y2, y3], label="Spring S2", linewidth=2, color=:green)
        plot!([x3, 0], [y3, 0], label="Spring S3", linewidth=2, color=:orange)
        xlims!(-2, 2)
        ylims!(-1, 4)
        title!("3-Point Mass-Spring System")
    end
    gif(anim, "mass_spring_system.gif", fps=50)
end

function main()
    se = Settings3Point()
    sol = simulate_system(se)
    animate_system(se, sol)
end

main()