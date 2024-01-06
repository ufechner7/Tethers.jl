# Tutorial example simulating a 3D mass-spring system with a nonlinear spring (no spring forces
# for l < l_0), n tether segments, reel-in and reel-out and continues callbacks. 
using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra, Timers, Parameters

@with_kw mutable struct Settings2 @deftype Float64
    g_earth::Vector{Float64} = [0.0, 0.0, -9.81] # gravitational acceleration     [m/s²]
    l0 = 50                                      # initial tether length             [m]
    v0 = 0                                       # initial speed                   [m/s]
    v_ro = 2                                     # reel-out speed                  [m/s]
    d_tether = 4                                 # tether diameter                  [mm]
    rho_tether = 724                             # density of Dyneema            [kg/m³]
    c_spring = 614600                            # unit spring constant              [N]
    damping = 473                                # unit damping constant            [Ns]
    segments::Int64 = 5                          # number of tether segments         [-]
    α0 = π/10                                    # initial tether angle            [rad]
    duration = 10                                # duration of the simulation        [s]
    dt = 0.02                                    # max time step in the solution     [s]
    save::Bool = false                           # save png files in folder video
    callbacks::Bool = true                       # use callbacks
end
                              
function calc_initial_state(se)
    POS0 = zeros(3, se.segments+1)
    VEL0 = zeros(3, se.segments+1)
    ACC0 = zeros(3, se.segments+1)
    SEGMENTS0 = zeros(3, se.segments) 
    UNIT_VECTORS0 = zeros(3, se.segments)
    for i in 1:se.segments+1
        l0 = -(i-1)*se.l0/se.segments
        v0 = (i-1)*se.v0/se.segments
        POS0[:, i] .= [sin(se.α0) * l0, 0, cos(se.α0) * l0]
        VEL0[:, i] .= [sin(se.α0) * v0, 0, cos(se.α0) * v0]
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
    @variables total_force(t)[1:3, 1:se.segments] = zeros(3, se.segments)
    D = Differential(t)

    eqs1 = vcat(D.(pos) ~ vel,
                D.(vel) ~ acc)
    eqs2 = []
    for i in se.segments:-1:1
        eqs2 = vcat(eqs2, segment[:, i] ~ pos[:, i+1] - pos[:, i])
        eqs2 = vcat(eqs2, norm1[i] ~ norm(segment[:, i]))
        eqs2 = vcat(eqs2, unit_vector[:, i] ~ -segment[:, i]/norm1[i])
        eqs2 = vcat(eqs2, rel_vel[:, i] ~ vel[:, i+1] - vel[:, i])
        eqs2 = vcat(eqs2, spring_vel[i] ~ -unit_vector[:, i] ⋅ rel_vel[:, i])
        eqs2 = vcat(eqs2, c_spr[i] ~ c_spring * (norm1[i] > length/se.segments))
        eqs2 = vcat(eqs2, spring_force[:, i] ~ (c_spr[i] * (norm1[i] - (length/se.segments)) + damping * spring_vel[i]) * unit_vector[:, i])
        if i == se.segments
            eqs2 = vcat(eqs2, total_force[:, i] ~ spring_force[:, i])
            eqs2 = vcat(eqs2, acc[:, i+1] .~ se.g_earth + total_force[:, i] / 0.5*(m_tether_particle))
        else
            eqs2 = vcat(eqs2, total_force[:, i] ~ spring_force[:, i]- spring_force[:, i+1])
            eqs2 = vcat(eqs2, acc[:, i+1] .~ se.g_earth + total_force[:, i] / m_tether_particle)
        end
    end
    eqs2 = vcat(eqs2, acc[:, 1] .~ zeros(3))
    eqs2 = vcat(eqs2, length ~ se.l0 + se.v_ro*t)
    eqs2 = vcat(eqs2, c_spring ~ se.c_spring / (length/se.segments))
    eqs2 = vcat(eqs2, m_tether_particle ~ mass_per_meter * (length/se.segments))
    eqs2 = vcat(eqs2, damping  ~ se.damping  / (length/se.segments))
    eqs = vcat(eqs1..., eqs2)

    if se.callbacks
        local cb
        for i in 1:se.segments
            cbi = [norm([pos[1, i+1] - pos[1, i], pos[2, i+1] - pos[2, i], pos[3, i+1] - pos[3, i]]) ~ abs(se.l0)/se.segments]
            if i == 1
                cb = cbi
            else
                cb = vcat(cb, cbi)
            end
        end
        @named sys = ODESystem(eqs, t; continuous_events = cb)
    else
        @named sys = ODESystem(eqs, t)
    end
    simple_sys = structural_simplify(sys)
    simple_sys, pos, vel
end

function simulate(se, simple_sys)
    tol = 1e-6
    tspan = (0.0, se.duration)
    ts    = 0:se.dt:se.duration
    prob = ODEProblem(simple_sys, nothing, tspan)
    @time sol = solve(prob, Rodas5(), dt=se.dt, abstol=tol, reltol=tol, saveat=ts)
    sol
end

function plot2d(se, sol, pos, reltime, line, sc, txt, j)
    index = findfirst(x -> x>=reltime, sol.t)
    x, z = Float64[], Float64[]
    for particle in 1:se.segments+1
        push!(x, (sol(sol.t, idxs=pos[1, particle]))[index])
        push!(z, (sol(sol.t, idxs=pos[3, particle]))[index])
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

function play(se, sol, pos)
    PyPlot.close()
    dt = 0.151
    ylim(-1.2*(se.l0+se.v_ro*se.duration), 0.5)
    xlim(-se.l0/2, se.l0/2)
    grid(true; color="grey", linestyle="dotted")
    tight_layout(rect=(0, 0, 0.98, 0.98))
    line, sc, txt = nothing, nothing, nothing
    start = time_ns()
    mkpath("video")
    for (j, time) in pairs(0:dt:se.duration)
        line, sc, txt = plot2d(se, sol, pos, time, line, sc, txt, j)
        wait_until(start + 0.5*time*1e9)
    end
    nothing
end

function main(se = Settings2();play_=true)
    simple_sys, pos, vel = model(se)
    sol = simulate(se, simple_sys)
    if play_
        play(se, sol, pos)
    end
    println("Events: $(Int64(round(length(sol.t)- se.duration/se.dt)-1))")
    sol, pos, vel
end

# plot vertical velocity of last particle with and without callbacks
function compare(se=Settings2())
    PyPlot.close()
    se.callbacks = false
    se.v_ro = 2
    se.v0 = 0
    sol, pos, vel = main(se, play_=false)
    vel_z1 = stack(sol[vel], dims=1)[:, 3, 5]
    plot(sol.t, vel_z1)
    se.callbacks = true
    sol, pos, vel = main(se, play_=false)
    vel_z2 = stack(sol[vel], dims=1)[:, 3, 5]
    plot(sol.t, vel_z2)
end

sol, pos, vel = main()
vel_z = stack(sol[vel], dims=1)[:, 3, 5]
nothing
