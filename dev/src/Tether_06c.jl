# Tutorial example simulating a 3D mass-spring system with a nonlinear spring (no spring forces
# for l < l_0), n tether segments, reel-in and reel-out and continues callbacks. 
using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra, Timers, Parameters
using ModelingToolkit: t_nounits as t, D_nounits as D
using ControlPlots

@with_kw mutable struct Settings2 @deftype Float64
    g_earth::Vector{Float64} = [0.0, 0.0, -9.81] # gravitational acceleration     [m/s²]
    l0 = 50                                      # initial tether length             [m]
    v0 = 0                                       # initial speed                   [m/s]
    v_ro = 2                                     # reel-out speed                  [m/s]
    d_tether = 4                                 # tether diameter                  [mm]
    rho_tether = 724                             # density of Dyneema            [kg/m³]
    c_spring = 614600                            # unit spring constant              [N]
    rel_compression_stiffness = 0.01             # relative compression stiffness    [-]    
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
    mass_per_meter = se.rho_tether * π * (se.d_tether/2000.0)^2
    @parameters c_spring0=se.c_spring/(se.l0/se.segments) l_seg=se.l0/se.segments
    @parameters rel_compression_stiffness = se.rel_compression_stiffness    
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
    @variables total_force(t)[1:3, 1:se.segments+1] = zeros(3, se.segments+1)

    eqs1 = vcat(D.(pos) .~ vel,
                D.(vel) .~ acc)
    eqs2 = vcat(eqs1...)

    for i in se.segments:-1:1
        eqs = [segment[:, i]        ~ pos[:, i+1] - pos[:, i],
               norm1[i]             ~ norm(segment[:, i]),
               unit_vector[:, i]    ~ -segment[:, i]/norm1[i],
               rel_vel[:, i]        ~ vel[:, i+1] - vel[:, i],
               spring_vel[i]        ~ -unit_vector[:, i] ⋅ rel_vel[:, i],
               c_spr[i]             ~ c_spring/(1+rel_compression_stiffness) 
                                      * (rel_compression_stiffness+(norm1[i] > length/se.segments)),               
               spring_force[:, i]  .~ (c_spr[i] * (norm1[i] - (length/se.segments)) 
                                       + damping * spring_vel[i]) * unit_vector[:, i]]
        if i == se.segments
            push!(eqs, total_force[:, i+1] .~ spring_force[:, i])
            push!(eqs, acc[:, i+1]         .~ se.g_earth + total_force[:, i+1] / (0.5*m_tether_particle))
            push!(eqs, total_force[:, i]   .~ spring_force[:, i-1] - spring_force[:, i])
        elseif i == 1
            push!(eqs, total_force[:, i]   .~ spring_force[:, i])
            push!(eqs, acc[:, i+1]         .~ se.g_earth + total_force[:, i+1] / m_tether_particle)            
        else
            push!(eqs, total_force[:, i]   .~ spring_force[:, i-1]- spring_force[:, i])
            push!(eqs, acc[:, i+1]         .~ se.g_earth + total_force[:, i+1] / m_tether_particle)
        end
        eqs2 = vcat(eqs2, reduce(vcat, eqs))
    end

    eqs = [acc[:, 1]         ~ zeros(3),
           length            ~ se.l0 + se.v_ro*t,
           c_spring          ~ se.c_spring / (length/se.segments),
           m_tether_particle ~ mass_per_meter * (length/se.segments),
           damping           ~ se.damping  / (length/se.segments)]
    eqs2 = vcat(eqs2, reduce(vcat, eqs))  

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
        @named sys = ODESystem(Symbolics.scalarize.(reduce(vcat, Symbolics.scalarize.(eqs2))), t; continuous_events = cb)
    else
        @named sys = ODESystem(Symbolics.scalarize.(reduce(vcat, Symbolics.scalarize.(eqs2))), t)
    end
    simple_sys = structural_simplify(sys)
    simple_sys, pos, vel
end

function simulate(se, simple_sys)
    tol = 1e-6
    tspan = (0.0, se.duration)
    ts    = 0:se.dt:se.duration
    prob = ODEProblem(simple_sys, nothing, tspan)
    elapsed_time = @elapsed sol = solve(prob, Rodas5P(), dt=se.dt, abstol=tol, reltol=tol, saveat=ts)
    elapsed_time = @elapsed sol = solve(prob, Rodas5P(), dt=se.dt, abstol=tol, reltol=tol, saveat=ts)
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
    i = 1; j = 0
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

function main(se = Settings2(); play_=true)
    simple_sys, pos, vel = model(se)
    sol, elapsed_time = simulate(se, simple_sys)
    println("Elapsed time: $(elapsed_time) s")
    if play_
        play(se, sol, pos)
    end
    println("Events: $(Int64(round(length(sol.t)- se.duration/se.dt)-1))")
    sol, pos, vel
end

# plot vertical velocity of last particle with and without callbacks
function compare(se=Settings2())
    se.callbacks = false
    se.v_ro = 2
    se.v0 = 0
    sol1, pos, vel = main(se, play_=false)
    vel_z1 = stack(sol1[vel], dims=1)[:, 3, 5]
    se.callbacks = true
    sol2, pos, vel = main(se, play_=false)
    vel_z2 = stack(sol2[vel], dims=1)[:, 3, 5]
    i, j = 1, 1
    vel1, vel2 = Float64[], Float64[]
    for time in sol1.t
        push!(vel1, vel_z1[i])
        while sol2.t[j] < time
            j += 1
        end
        push!(vel2, vel_z2[j])
        i += 1
    end
    plot(sol1.t, [vel1, vel2]; xlabel="time [s]", ylabel="velocity [m/s]", labels=["without callbacks" "with callbacks"]) 
end

if (! @isdefined __BENCH__) || __BENCH__ == false
    sol, pos, vel = main()
    vel_z = stack(sol[vel], dims=1)[:, 3, 5]    
end
__BENCH__ = false
nothing
