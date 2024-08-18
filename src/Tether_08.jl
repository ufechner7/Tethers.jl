# Tutorial example simulating a 3D mass-spring system with a nonlinear spring (1% stiffnes
# for l < l_0), n tether segments, tether drag and reel-in and reel-out. 
# New feature: A steady state solver shall be used to allow different initial conditions.
using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra, Timers, Parameters
using ModelingToolkit: t_nounits as t, D_nounits as D
using ControlPlots

@with_kw mutable struct Settings3 @deftype Float64
    g_earth::Vector{Float64} = [0.0, 0.0, -9.81] # gravitational acceleration     [m/s²]
    v_wind_tether::Vector{Float64} = [2, 0.0, 0.0]
    rho = 1.225
    cd_tether = 0.958
    l0 = 50                                      # initial tether length             [m]
    v_ro = 2                                     # reel-out speed                  [m/s]
    d_tether = 4                                 # tether diameter                  [mm]
    rho_tether = 724                             # density of Dyneema            [kg/m³]
    c_spring = 614600                            # unit spring constant              [N]
    rel_compression_stiffness = 0.01             # relative compression stiffness    [-]
    damping = 473                                # unit damping constant            [Ns]
    segments::Int64 = 5                          # number of tether segments         [-]
    α0 = π/10                                    # initial tether angle            [rad]
    duration = 10                                # duration of the simulation        [s]
    save::Bool = false                           # save png files in folder video
end

function set_tether_diameter!(se, d; c_spring_4mm = 614600, damping_4mm = 473)
    se.d_tether = d
    se.c_spring = c_spring_4mm * (d/4.0)^2
    se.damping = damping_4mm * (d/4.0)^2
end
                              
function calc_initial_state(se; p1, p2)
    # start with a linear interpolation between p1 and p2
    if isnothing(p2)
        z  = cos(se.α0) * se.l0
        y  = sin(se.α0) * se.l0
        p2 = [p1[1], p1[2] - y, p1[3] - z]
        println("p2: ", p2)
    end
    POS0 = zeros(3, se.segments+1)
    VEL0 = zeros(3, se.segments+1)
    ACC0 = zeros(3, se.segments+1)
    SEGMENTS0 = zeros(3, se.segments) 
    UNIT_VECTORS0 = zeros(3, se.segments)
    for i in 1:se.segments+1
        Δ = (p2-p1) / se.segments
        POS0[:, i] .= p1 + (i-1) * Δ
    end
    for i in 1:se.segments
        ACC0[:, i+1] .= se.g_earth
        UNIT_VECTORS0[:, i] .= [0, 0, 1.0]
        SEGMENTS0[:, i] .= POS0[:, i+1] - POS0[:, i]
    end
    POS0, VEL0, ACC0, SEGMENTS0, UNIT_VECTORS0
end

function model(se; p1=[0,0,0], p2=nothing, fix_p1=true, fix_p2=false)
    if !isnothing(p1)
        @assert isa(p1, AbstractVector) || error("p1 must be a vector")
        @assert (length(p1) == 3)       || error("p1 must have length 3")
    end
    if !isnothing(p2)
        @assert isa(p2, AbstractVector) || error("p2 must be a vector")
        @assert (length(p2) == 3)       || error("p2 must have length 3")
    end
    POS0, VEL0, ACC0, SEGMENTS0, UNIT_VECTORS0 = calc_initial_state(se; p1, p2)
    mass_per_meter = se.rho_tether * π * (se.d_tether/2000.0)^2
    @parameters c_spring0=se.c_spring/(se.l0/se.segments) l_seg=se.l0/se.segments
    @parameters rel_compression_stiffness = se.rel_compression_stiffness
    @variables pos(t)[1:3, 1:se.segments+1]  = POS0
    @variables vel(t)[1:3, 1:se.segments+1]  = VEL0
    @variables acc(t)[1:3, 1:se.segments+1]  = ACC0
    @variables segment(t)[1:3, 1:se.segments]  = SEGMENTS0
    @variables unit_vector(t)[1:3, 1:se.segments]  = UNIT_VECTORS0
    @variables len(t) = se.l0
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
    @variables total_force(t)[1:3, 1:se.segments+1] = zeros(3, se.segments+1)

    # basic differential equations
    eqs1 = vcat(D.(pos) .~ vel,
                D.(vel) .~ acc)
    eqs2 = vcat(eqs1...)
    # loop over all segments to calculate the spring and drag forces
    for i in 1:se.segments
        eqs = [segment[:, i]      ~ pos[:, i+1] - pos[:, i],
               norm1[i]           ~ norm(segment[:, i]),
               unit_vector[:, i]  ~ -segment[:, i]/norm1[i],
               rel_vel[:, i]      ~ vel[:, i+1] - vel[:, i],
               spring_vel[i]      ~ -unit_vector[:, i] ⋅ rel_vel[:, i],
               c_spr[i]           ~ c_spring/(1+rel_compression_stiffness) 
                                     * (rel_compression_stiffness+(norm1[i] > len/se.segments)),
               spring_force[:, i] ~ (c_spr[i] * (norm1[i] - (len/se.segments)) 
                                     + damping * spring_vel[i]) * unit_vector[:, i],
               v_apparent[:, i]   ~ se.v_wind_tether .- (vel[:, i] + vel[:, i+1])/2,
               v_app_perp[:, i]   ~ v_apparent[:, i] - (v_apparent[:, i] ⋅ unit_vector[:, i]) .* unit_vector[:, i],
               norm_v_app[i]      ~ norm(v_app_perp[:, i]),
               half_drag_force[:, i] ~ 0.25 * se.rho * se.cd_tether * norm_v_app[i] * (norm1[i]*se.d_tether/1000.0)
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
    eqs = [len            ~ se.l0 + se.v_ro*t,
           c_spring          ~ se.c_spring / (len/se.segments),
           m_tether_particle ~ mass_per_meter * (len/se.segments),
           damping           ~ se.damping  / (len/se.segments)]
    eqs2 = vcat(eqs2, reduce(vcat, eqs))  
        
    @named sys = ODESystem(Symbolics.scalarize.(reduce(vcat, Symbolics.scalarize.(eqs2))), t)
    simple_sys = structural_simplify(sys)
    simple_sys, pos, vel, len, c_spr
end

function simulate(se, simple_sys)
    dt = 0.02
    tol = 1e-6
    tspan = (0.0, se.duration)
    ts    = 0:dt:se.duration
    prob = ODEProblem(simple_sys, nothing, tspan)
    elapsed_time = @elapsed sol = solve(prob, KenCarp4(autodiff=false); dt, abstol=tol, reltol=tol, saveat=ts)
    elapsed_time = @elapsed sol = solve(prob, KenCarp4(autodiff=false); dt, abstol=tol, reltol=tol, saveat=ts)
    sol, elapsed_time
end

function play(se, sol, pos)
    dt = 0.151
    ylim = (-1.2 * (se.l0 + se.v_ro*se.duration), 0.5)
    xlim = (-se.l0/2, se.l0/2)
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
        plot2d(sol[pos][i], time; segments=se.segments, xlim, ylim, xy)
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
    play(se, sol, pos)
    println("Elapsed time: $(elapsed_time) s, speed: $(round(se.duration/elapsed_time)) times real-time")
    println("Number of evaluations per step: ", round(sol.stats.nf/(se.duration/0.02), digits=1))
    sol, pos, vel, simple_sys
end

# if (! @isdefined __BENCH__) || __BENCH__ == false
#     sol, pos, vel, simple_sys = main()
# end
# __BENCH__ = false
nothing
