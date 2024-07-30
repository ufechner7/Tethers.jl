# Tutorial example simulating a 3D mass-spring system with a nonlinear spring (no spring forces
# for l < l_0), n tether segments and reel-in and reel-out. 
using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra, Timers, Parameters

# TODO: Add aerodynamic drag

@with_kw mutable struct Settings @deftype Float64
    g_earth::Vector{Float64} = [0.0, 0.0, -9.81] # gravitational acceleration     [m/s²]
    v_wind_tether::Vector{Float64} = [5.0, 0.0, 0.0]
    rho = 1.225
    cd_tether = 0.958
    l0 = 50                                      # initial tether length             [m]
    v_ro = 2                                     # reel-out speed                  [m/s]
    d_tether = 10                                 # tether diameter                  [mm]
    rho_tether = 724                             # density of Dyneema            [kg/m³]
    c_spring = 614600                            # unit spring constant              [N]
    damping = 473                                # unit damping constant            [Ns]
    segments::Int64 = 5                          # number of tether segments         [-]
    α0 = π/10                                    # initial tether angle            [rad]
    duration = 10.0                             # duration of the simulation        [s]
    save::Bool = false                           # save png files in folder video
    dt = 0.2
end
                              
function calc_initial_state(se)
    POS0 = zeros(3, se.segments+1)
    VEL0 = zeros(3, se.segments+1)
    ACC0 = zeros(3, se.segments+1)
    SEGMENTS0 = zeros(3, se.segments) 
    UNIT_VECTORS0 = zeros(3, se.segments)
    for i in 1:se.segments+1
        l0 = -(i-1)*se.l0/se.segments
        POS0[:, i] .= [sin(se.α0) * l0, 0, cos(se.α0) * l0]
        VEL0[:, i] .= [0, 0, 0]
    end
    for i in 1:se.segments
        ACC0[:, i+1] .= se.g_earth
        UNIT_VECTORS0[:, i] .= [0, 0, 1.0]
        SEGMENTS0[:, i] .= POS0[:, i+1] - POS0[:, i]
    end
    POS0, VEL0, ACC0, SEGMENTS0, UNIT_VECTORS0
end

function model(se)
    POS0, VEL0, ACC0, SEGMENTS0, UNIT_VECTORS0 = calc_initial_state(se)
    mass_per_meter = se.rho_tether * se.segments * (se.d_tether/2000.0)^2
    @parameters c_spring0=se.c_spring/(se.l0/se.segments) l_seg=se.l0/se.segments
    @variables t 
    @variables pos(t)[1:3, 1:se.segments+1]  = POS0
    @variables vel(t)[1:3, 1:se.segments+1]  = VEL0
    @variables acc(t)[1:3, 1:se.segments+1]  = ACC0
    @variables segment(t)[1:3, 1:se.segments]  = SEGMENTS0
    @variables unit_vector(t)[1:3, 1:se.segments]  = UNIT_VECTORS0
    @variables length(t) = se.l0
    @variables c_spring(t) = c_spring0
    @variables damping(t) = se.damping  / l_seg
    @variables m_tether_particle(t) = mass_per_meter * l_seg
    @variables norm1(t)[1:se.segments] = l_seg * ones(se.segments)
    @variables rel_vel(t)[1:3, 1:se.segments]  = zeros(3, se.segments)
    @variables spring_vel(t)[1:se.segments] = zeros(se.segments)
    @variables c_spr(t)[1:se.segments] = c_spring0 * ones(se.segments)
    @variables spring_force(t)[1:3, 1:se.segments] = zeros(3, se.segments)
    @variables v_apparent(t)[1:3, 1:se.segments] = zeros(3, se.segments)
    @variables v_app_perp(t)[1:3, 1:se.segments] = zeros(3, se.segments)
    @variables norm_v_app(t)[1:se.segments] = ones(se.segments)
    @variables half_drag_force(t)[1:3, 1:se.segments] = zeros(3, se.segments)
    @variables total_force(t)[1:3, 1:se.segments] = zeros(3, se.segments)

    D = Differential(t)

    eqs1 = vcat(D.(pos) .~ vel,
                D.(vel) .~ acc)
    eqs2 = []
    for i in se.segments:-1:1
        eqs2 = vcat(eqs2, segment[:, i] .~ pos[:, i+1] - pos[:, i])
        eqs2 = vcat(eqs2, norm1[i] ~ norm(segment[:, i]))
        eqs2 = vcat(eqs2, unit_vector[:, i] .~ -segment[:, i]/norm1[i])
        eqs2 = vcat(eqs2, rel_vel[:, i] .~ vel[:, i+1] - vel[:, i])
        eqs2 = vcat(eqs2, spring_vel[i] .~ -unit_vector[:, i] ⋅ rel_vel[:, i])
        eqs2 = vcat(eqs2, c_spr[i] .~ c_spring * (norm1[i] > length/se.segments))
        eqs2 = vcat(eqs2, spring_force[:, i] .~ (c_spr[i] * (norm1[i] - (length/se.segments)) + damping * spring_vel[i]) * unit_vector[:, i])

        eqs2 = vcat(eqs2, v_apparent[:, i] .~ se.v_wind_tether .- (vel[:, i] + vel[:, i+1])/2)
        eqs2 = vcat(eqs2, v_app_perp[:, i] .~ v_apparent[:, i] - (v_apparent[:, i] ⋅ unit_vector[:, i]) .* unit_vector[:, i])
        eqs2 = vcat(eqs2, norm_v_app[i] ~ norm(v_app_perp[:, i]))
        eqs2 = vcat(eqs2, half_drag_force[:, i] .~ (0.25 * se.rho * se.cd_tether * norm_v_app[i] * (norm1[i]*se.d_tether/1000.0)) .* v_app_perp[:, i])
        if i == se.segments
            eqs2 = vcat(eqs2, total_force[:, i] .~ spring_force[:, i] + half_drag_force[:,i] + half_drag_force[:,i-1])
            eqs2 = vcat(eqs2, acc[:, i+1] .~ se.g_earth + total_force[:, i] / 0.5*(m_tether_particle))
        elseif i == 1
            eqs2 = vcat(eqs2, total_force[:, i] .~ spring_force[:, i]- spring_force[:, i+1] + half_drag_force[:,i])
            eqs2 = vcat(eqs2, acc[:, i+1] .~ se.g_earth + total_force[:, i] / m_tether_particle)
        else
            # println(spring_force[:, i]+drag_force[:, i])
            # println(drag_force[:, i])
            eqs2 = vcat(eqs2, total_force[:, i] .~ spring_force[:, i]- spring_force[:, i+1] + half_drag_force[:,i] + half_drag_force[:,i-1])
            eqs2 = vcat(eqs2, acc[:, i+1] .~ se.g_earth + total_force[:, i] / m_tether_particle)
        end
    end
    eqs2 = vcat(eqs2, acc[:, 1] .~ zeros(3))
    eqs2 = vcat(eqs2, length .~ se.l0 + se.v_ro*t)
    eqs2 = vcat(eqs2, c_spring .~ se.c_spring / (length/se.segments))
    eqs2 = vcat(eqs2, m_tether_particle .~ mass_per_meter * (length/se.segments))
    eqs2 = vcat(eqs2, damping  .~ se.damping  / (length/se.segments))
    eqs = vcat(eqs1..., eqs2)
        
    @named sys = ODESystem(eqs, t)
    simple_sys = structural_simplify(sys)
    simple_sys, pos, vel
end

function next_step(se, integrator)
    @time step!(integrator, se.dt, true)
end

function plot2d(se, integrator, pos, reltime, line, sc, txt, j)
    index = Int64(round(reltime*50+1))
    x, z = Float64[], Float64[]
    for particle in 1:se.segments+1
        push!(x, integrator.u[1:6:end][particle])
        push!(z, integrator.u[3:6:end][particle])
    end
    z_max = maximum(z)
    if isnothing(line)
        line, = plot(x,z; linewidth="1")
        sc  = scatter(x, z; s=15, color="red") 
        txt = annotate("t=$(round(reltime,digits=1)) s",  
                        xy=(se.l0/4.2, z_max-7), fontsize = 12)
    else
        line.set_xdata(x)
        line.set_ydata(z)
        sc.set_offsets(hcat(x,z))
        txt.set_text("t=$(round(reltime,digits=1)) s")
        gcf().canvas.draw()
    end
    if se.save
        PyPlot.savefig("video/"*"img-"*lpad(j,4,"0"))
    end
    line, sc, txt
end

function play(se, simple_sys, pos)
    PyPlot.close()
    ylim(-1.2*(se.l0+se.v_ro*se.duration), 0.5)
    xlim(-se.l0/2, se.l0/2)
    grid(true; color="grey", linestyle="dotted")
    tight_layout(rect=(0, 0, 0.98, 0.98))
    line, sc, txt = nothing, nothing, nothing
    start = time_ns()
    mkpath("video")
    tspan = (0.0, se.dt)
    prob = ODEProblem(simple_sys, nothing, tspan)
    tol=1e-6
    integrator = init(prob, Rodas5P(); dt=se.dt, abstol=tol, reltol=tol, save_everystep=false)
    for (j, time) in pairs(0:se.dt:se.duration)
        next_step(se, integrator)
        println(integrator.u[4:6:end])
        line, sc, txt = plot2d(se, integrator, pos, time, line, sc, txt, j)
        wait_until(start + 1.0*time*1e9)
    end
    nothing
end

function main()
    se = Settings()
    simple_sys, pos, vel = model(se)
    # println("sol and pos ", sol, "\n\t", pos)
    play(se, simple_sys, pos)
end

main()
