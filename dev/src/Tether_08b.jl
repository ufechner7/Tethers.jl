# Tutorial example simulating a 3D mass-spring system with a nonlinear spring (1% stiffnes
# for l < l_0), n tether segments, tether drag and reel-in and reel-out. 
# New feature: A steady state solver is used to find the initial tether shape for any
# given pair of endpoints, which is then used as the initial condition for the simulation.
using ModelingToolkit, OrdinaryDiffEq, SteadyStateDiffEq, LinearAlgebra, Timers, Parameters, ControlPlots
tic()
using ModelingToolkit: t_nounits as t, D_nounits as D
using ControlPlots

@with_kw mutable struct Settings3 @deftype Float64
    g_earth::Vector{Float64} = [0.0, 0.0, -9.81] # gravitational acceleration     [m/s²]
    v_wind_tether::Vector{Float64} = [4, 0.0, 0.0]
    rho = 1.225
    cd_tether = 0.958
    l0 = 70                                      # initial tether length             [m]
    v_ro = 1.2                                   # reel-out speed                  [m/s]
    d_tether = 12                                 # tether diameter                  [mm]
    rho_tether = 724                             # density of Dyneema            [kg/m³]
    c_spring = 614600                            # unit spring constant              [N]
    rel_compression_stiffness = 0.01             # relative compression stiffness    [-]
    damping = 473                                # unit damping constant            [Ns]
    segments::Int64 = 10                         # number of tether segments         [-]
    α0 = π/10                                    # initial tether angle            [rad]
    duration = 30                                # duration of the simulation        [s]
    save::Bool = false                           # save png files in folder video
end

function set_tether_diameter!(se, d; c_spring_4mm = 614600, damping_4mm = 473)
    se.d_tether = d
    se.c_spring = c_spring_4mm * (d/4.0)^2
    se.damping = damping_4mm * (d/4.0)^2
end
                              
function calc_initial_state(se; p1, p2)
    # calculate p2 based on se.α0 and se.l0 if not given
    if isnothing(p2)
        z  = cos(se.α0) * se.l0
        y  = sin(se.α0) * se.l0
        p2 = [p1[1], p1[2] - y, p1[3] - z]
        println("p2: ", p2)
    end
    POS0 = zeros(3, se.segments+1)
    VEL0 = zeros(3, se.segments+1)
    # use a linear interpolation between p1 and p2 for the intermediate points
    for i in 1:se.segments+1
        Δ = (p2-p1) / se.segments
        POS0[:, i] .= p1 + (i-1) * Δ
    end
    POS0, VEL0
end

function model(se; p1=[0,0,0], p2=nothing, fix_p1=true, fix_p2=false)
    if ! isnothing(p1)
        @assert isa(p1, AbstractVector) || error("p1 must be a vector")
        @assert (length(p1) == 3)       || error("p1 must have length 3")
    else
        @assert ! fix_p1                || error("if p1 undefined it cannot be fixed")
    end

    if ! isnothing(p2)
        @assert isa(p2, AbstractVector) || error("p2 must be a vector")
        @assert (length(p2) == 3)       || error("p2 must have length 3")
    else
        @assert ! fix_p2                || error("if p2 undefined it cannot be fixed")
    end
    # straight line approximation for the tether
    POS0, VEL0 = calc_initial_state(se; p1, p2)
    # find steady state
    v_ro = se.v_ro      # save the reel-out speed
    se.v_ro = 0         # v_ro must be zero, otherwise finding the steady state is not possible
    simple_sys, pos, =  model(se, p1, p2, true, true, POS0, VEL0)
    tspan = (0.0, se.duration)
    prob = ODEProblem(simple_sys, nothing, tspan)
    prob1 = SteadyStateProblem(prob)
    sol1 = solve(prob1, DynamicSS(KenCarp4(autodiff=false)))
    POS0 = sol1[pos]
    # create the real model
    se.v_ro = v_ro
    model(se, p1, p2, fix_p1, fix_p2, POS0, VEL0)
end
function model(se, p1, p2, fix_p1, fix_p2, POS0, VEL0)
    mass_per_meter = se.rho_tether * π * (se.d_tether/2000.0)^2
    @parameters c_spring0=se.c_spring/(se.l0/se.segments) l_seg=se.l0/se.segments
    @parameters rel_compression_stiffness = se.rel_compression_stiffness
    @variables begin 
        pos(t)[1:3, 1:se.segments+1]  = POS0
        vel(t)[1:3, 1:se.segments+1]  = VEL0
        acc(t)[1:3, 1:se.segments+1]
        segment(t)[1:3, 1:se.segments]
        unit_vector(t)[1:3, 1:se.segments]
        l_spring(t), c_spring(t), damping(t), m_tether_particle(t)
        len(t)[1:se.segments]
        rel_vel(t)[1:3, 1:se.segments]
        spring_vel(t)[1:se.segments]
        c_spr(t)[1:se.segments]
        spring_force(t)[1:3, 1:se.segments]
        v_apparent(t)[1:3, 1:se.segments]
        v_app_perp(t)[1:3, 1:se.segments]
        norm_v_app(t)[1:se.segments]
        half_drag_force(t)[1:3, 1:se.segments]
        total_force(t)[1:3, 1:se.segments+1]
    end

    # basic differential equations
    eqs1 = vcat(D.(pos) .~ vel,
                D.(vel) .~ acc)
    eqs2 = vcat(eqs1...)
    # loop over all segments to calculate the spring and drag forces
    for i in 1:se.segments
        eqs = [segment[:, i]      ~ pos[:, i+1] - pos[:, i],
               len[i]             ~ norm(segment[:, i]),
               unit_vector[:, i]  ~ -segment[:, i]/len[i],
               rel_vel[:, i]      ~ vel[:, i+1] - vel[:, i],
               spring_vel[i]      ~ -unit_vector[:, i] ⋅ rel_vel[:, i],
               c_spr[i]           ~ c_spring / (1+rel_compression_stiffness) 
                                     * (rel_compression_stiffness+(len[i] > l_spring)),
               spring_force[:, i] ~ (c_spr[i] * (len[i] - l_spring) 
                                     + damping * spring_vel[i]) * unit_vector[:, i],
               v_apparent[:, i]   ~ se.v_wind_tether .- (vel[:, i] + vel[:, i+1])/2,
               v_app_perp[:, i]   ~ v_apparent[:, i] - (v_apparent[:, i] ⋅ unit_vector[:, i]) .* unit_vector[:, i],
               norm_v_app[i]      ~ norm(v_app_perp[:, i]),
               half_drag_force[:, i] ~ 0.25 * se.rho * se.cd_tether * norm_v_app[i] * (len[i]*se.d_tether/1000.0)
                                        * v_app_perp[:, i]]
        eqs2 = vcat(eqs2, reduce(vcat, eqs))
    end
    # loop over all tether particles to apply the forces and calculate the accelerations
    for i in 1:(se.segments+1)
        eqs = []
        if i == se.segments+1
            push!(eqs, total_force[:, i] ~ spring_force[:, i-1] + half_drag_force[:, i-1])
            if isnothing(p2) || ! fix_p2
                push!(eqs, acc[:, i]         ~ se.g_earth + total_force[:, i] / (0.5 * m_tether_particle))
            else
                push!(eqs, acc[:, i]         ~ zeros(3))
            end
        elseif i == 1
            push!(eqs, total_force[:, i] ~ spring_force[:, i] + half_drag_force[:, i])
            if isnothing(p1) || ! fix_p1
                push!(eqs, acc[:, i]     ~ se.g_earth + total_force[:, i] / (0.5 * m_tether_particle))
            else
                push!(eqs, acc[:, i]     ~ zeros(3))
            end
        else
            push!(eqs, total_force[:, i] ~ spring_force[:, i-1] - spring_force[:, i] 
                                           + half_drag_force[:, i-1] + half_drag_force[:, i])
            push!(eqs, acc[:, i]         ~ se.g_earth + total_force[:, i] / m_tether_particle)
        end
        eqs2 = vcat(eqs2, reduce(vcat, eqs))
    end
    # scalar equations
    eqs = [l_spring   ~ (se.l0 + se.v_ro*t) / se.segments,
    c_spring          ~ se.c_spring / l_spring,
    m_tether_particle ~ mass_per_meter * l_spring,
    damping           ~ se.damping  / l_spring]
    eqs2 = vcat(eqs2, reduce(vcat, eqs))  
        
    @named sys = ODESystem(reduce(vcat, Symbolics.scalarize.(eqs2)), t)
    simple_sys = structural_simplify(sys)
    simple_sys, pos, vel, len, c_spr
end

function simulate(se, simple_sys)
    dt = 0.02
    tol = 1e-6
    tspan = (0.0, se.duration)
    ts    = 0:dt:se.duration
    prob = ODEProblem(simple_sys, nothing, tspan)
    toc()
    elapsed_time = @elapsed sol = solve(prob, FBDF(autodiff=true); dt, abstol=tol, reltol=tol, saveat=ts)
    elapsed_time = @elapsed sol = solve(prob, FBDF(autodiff=true); dt, abstol=tol, reltol=tol, saveat=ts)
    sol, elapsed_time
end

function play(se, sol, pos)
    dt = 0.05
    ylim = (-1.2 * (se.l0 + se.v_ro*se.duration), 0.5)
    xlim = (-se.l0, se.l0)
    mkpath("video")
    z_max = 0.0
    # text position
    xy = (se.l0/4.2, z_max-7)
    start = time_ns()
    i = 1; j = 0
    for time in 0:dt:se.duration
        # while we run the simulation in steps of 20ms, we update the plot only every 150ms
        # therefore we have to skip some steps of the result
        while sol.t[i] < time
            i += 1
        end
        plot2d(sol[pos][i], time; segments=se.segments, xlim, ylim, xy, fig="Tether_08", figsize=(12.8, 9.6))
        if se.save
            ControlPlots.plt.savefig("video/"*"img-"*lpad(j, 4, "0"))
        end
        j += 1
        wait_until(start + 0.5 * time * 1e9)
    end
    if se.save
        include("export_gif.jl")
    end
    nothing
end

function main(; p1=[0,0,0], p2=nothing, fix_p1=true, fix_p2=false)
    global sol, pos, vel, len, c_spr
    se = Settings3()
    set_tether_diameter!(se, se.d_tether) # adapt spring and damping constants to tether diameter
    simple_sys, pos, vel, len, c_spr = model(se; p1, p2, fix_p1, fix_p2)
    sol, elapsed_time = simulate(se, simple_sys)
    if @isdefined __PC
        return sol, pos, vel, simple_sys
    end
    play(se, sol, pos)
    println("Elapsed time: $(elapsed_time) s, speed: $(round(se.duration/elapsed_time)) times real-time")
    println("Number of evaluations per step: ", round(sol.stats.nf/(se.duration/0.02), digits=1))
    sol, pos, vel, simple_sys
end

main(p2=[-40,0,-47], fix_p2=false);

nothing

