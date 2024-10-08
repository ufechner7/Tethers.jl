# Tutorial example simulating a 3D mass-spring system with a nonlinear spring (no spring forces
# for l < l_0), n tether segments and reel-in and reel-out. Can create a video of the simulation.
# For creating the video, set save=true in the Settings struct.
using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra, Timers, Parameters, ControlPlots
using ModelingToolkit: t_nounits as t, D_nounits as D
using ControlPlots

@with_kw mutable struct Settings @deftype Float64
    g_earth::Vector{Float64} = [0.0, 0.0, -9.81] # gravitational acceleration     [m/s²]
    l0 = 50.0                                    # initial tether length             [m]
    v_ro = 2                                     # reel-out speed                  [m/s]
    d_tether = 4                                 # tether diameter                  [mm]
    rho_tether = 724                             # density of Dyneema            [kg/m³]
    c_spring = 614600                            # unit spring constant              [N]
    damping = 473                                # unit damping constant            [Ns]
    segments::Int64 = 5                          # number of tether segments         [-]
    α0 = π/10                                    # initial tether angle            [rad]
    duration = 10                                # duration of the simulation        [s]
    save::Bool = false                           # save png files in folder video
end
                              
function calc_initial_state(se)
    POS0 = zeros(3, se.segments+1)
    VEL0 = zeros(3, se.segments+1)
    for i in 1:se.segments+1
        l0 = -(i-1)*se.l0/se.segments
        POS0[:, i] .= [sin(se.α0) * l0, 0, cos(se.α0) * l0]
        VEL0[:, i] .= [0, 0, 0]
    end
    POS0, VEL0
end

function model(se)
    POS0, VEL0 = calc_initial_state(se)
    mass_per_meter = se.rho_tether * π * (se.d_tether/2000.0)^2
    @parameters c_spring0=se.c_spring/(se.l0/se.segments) l_seg=se.l0/se.segments
    @variables pos(t)[1:3, 1:se.segments+1]  = POS0
    @variables vel(t)[1:3, 1:se.segments+1]  = VEL0
    @variables acc(t)[1:3, 1:se.segments+1]
    @variables segment(t)[1:3, 1:se.segments]
    @variables unit_vector(t)[1:3, 1:se.segments]
    @variables len(t)
    @variables c_spring(t)
    @variables damping(t)
    @variables m_tether_particle(t)
    @variables norm1(t)[1:se.segments]
    @variables rel_vel(t)[1:3, 1:se.segments]
    @variables spring_vel(t)[1:se.segments]
    @variables c_spr(t)[1:se.segments]
    @variables spring_force(t)[1:3, 1:se.segments]
    @variables total_force(t)[1:3, 1:se.segments+1]

    # basic differential equations
    eqs1 = vcat(D.(pos) .~ vel,
                D.(vel) .~ acc)
    eqs2 = vcat(eqs1...)

    # loop over all segments to calculate the spring forces
    for i in 1:se.segments
        eqs = [segment[:, i]      ~ pos[:, i+1] - pos[:, i],
               norm1[i]           ~ norm(segment[:, i]),
               unit_vector[:, i]  ~ -segment[:, i]/norm1[i],
               rel_vel[:, i]      ~ vel[:, i+1] - vel[:, i],
               spring_vel[i]      ~ -unit_vector[:, i] ⋅ rel_vel[:, i],
               c_spr[i]           ~ c_spring/1.01 * (0.01 + (norm1[i] > len/se.segments)),
               spring_force[:, i] ~ (c_spr[i] * (norm1[i] - (len/se.segments)) + damping * spring_vel[i]) * unit_vector[:, i]]
        eqs2 = vcat(eqs2, reduce(vcat, eqs))
    end
    # loop over all tether particles to apply the forces and calculate the accelerations
    for i in 1:(se.segments+1)
        eqs = []
        if i == se.segments+1
            push!(eqs, total_force[:, i] ~ spring_force[:, i-1])
            push!(eqs, acc[:, i]         ~ se.g_earth + total_force[:, i] / (0.5 * m_tether_particle))
        elseif i == 1
            push!(eqs, total_force[:, i] ~ spring_force[:, i])
            push!(eqs, acc[:, i]         ~ zeros(3))
        else
            push!(eqs, total_force[:, i] ~ spring_force[:, i-1] - spring_force[:, i] )
            push!(eqs, acc[:, i]         ~ se.g_earth + total_force[:, i] / m_tether_particle)
        end
        eqs2 = vcat(eqs2, reduce(vcat, eqs))
    end
    # scalar equations
    eqs = [len               ~ se.l0 + se.v_ro*t,
           c_spring          ~ se.c_spring / (len/se.segments),
           m_tether_particle ~ mass_per_meter * (len/se.segments),
           damping           ~ se.damping  / (len/se.segments)]
    eqs2 = vcat(eqs2, reduce(vcat, eqs))  
        
    @named sys = ODESystem(reduce(vcat, Symbolics.scalarize.(eqs2)), t)
    simple_sys = structural_simplify(sys)
    simple_sys, pos, vel
end

function simulate(se, simple_sys)
    dt = 0.02
    tol = 1e-6
    tspan = (0.0, se.duration)
    ts    = 0:dt:se.duration
    prob = ODEProblem(simple_sys, nothing, tspan)
    elapsed_time = @elapsed sol = solve(prob, FBDF(), dt=dt, abstol=tol, reltol=tol, saveat=ts)
    sol, elapsed_time
end

function play(se, sol, pos)
    dt = 0.151
    ylim = (-1.2*(se.l0+se.v_ro*se.duration), 0.5)
    xlim = (-se.l0/2, se.l0/2)
    mkpath("video")
    z_max = 0.0
    # text position
    xy = (se.l0/4.2, z_max-7)
    start = time_ns()
    i = 1
    j = 0
    for time in 0:dt:se.duration
        # while we run the simulation in steps of 20ms, we update the plot only every 150ms
        # therefore we have to skip some steps of the result
        while sol.t[i] < time
            i += 1
        end
        plot2d(sol[pos][i], time; segments=se.segments, xlim, ylim, xy)
        if se.save
            ControlPlots.plt.savefig("video/"*"img-"*lpad(j,4,"0"))
        end
        j += 1
        wait_until(start + 0.5*time*1e9)
    end
    if se.save
        include("export_gif.jl")
    end
    nothing
end

function main()
    se = Settings()
    simple_sys, pos, vel = model(se)
    sol, elapsed_time = simulate(se, simple_sys)
    play(se, sol, pos)
    println("Elapsed time: $(elapsed_time) s, speed: $(round(se.duration/elapsed_time)) times real-time")
    println("Number of evaluations per step: ", round(sol.stats.nf/(se.duration/0.02), digits=1))
end

main()
