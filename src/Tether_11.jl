# Tutorial example simulating a 3D mass-spring system with a nonlinear spring (1% stiffnes
# for l < l_0), n tether segments, tether drag and reel-in and reel-out. 
# New feature: A steady state solver is used to find the initial tether shape for any
# given pair of endpoints, which is then used as the initial condition for the simulation.
using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra, Timers, Parameters, ControlPlots, Optimization, OptimizationOptimJL
tic()
using ModelingToolkit: t_nounits as t, D_nounits as D
using ControlPlots

@with_kw mutable struct Settings3 @deftype Float64
    g_earth::Vector{Float64} = [0.0, 0.0, -9.81] # gravitational acceleration     [m/s²]
    v_wind_tether::Vector{Float64} = [10.0, 0.0, 0.0]
    rho = 1.225
    cd_tether = 0.958
    l0 = 70                                      # initial tether length             [m]
    v_ro = 0.3                                   # reel-out speed                  [m/s]
    d_tether = 4                                 # tether diameter                  [mm]
    rho_tether = 724                             # density of Dyneema            [kg/m³]
    c_spring = 614600                            # unit spring constant              [N]
    rel_compression_stiffness = 0.01             # relative compression stiffness    [-]
    damping = 473                                # unit damping constant            [Ns]
    segments::Int64 = 2                         # number of tether segments         [-]
    α0 = π/10                                    # initial tether angle            [rad]
    duration = 30                                # duration of the simulation        [s]
    save::Bool = false                           # save png files in folder video
end

function set_tether_diameter!(se, d; c_spring_4mm = 614600, damping_4mm = 473)
    se.d_tether = d
    se.c_spring = c_spring_4mm * (d/4.0)^2
    se.damping = damping_4mm * (d/4.0)^2
end

"""
Use an optimizationproblem to find pos0 and vel0
"""
function calc_initial_state(se, ss, eqs, p1, elevation, tether_length)
    global prob
    # remove the differential equations
    diff_eqs = equations(ss)
    rm_idxs = []
    for i in eachindex(eqs)
        for d in diff_eqs
            if isequal(d.lhs, eqs[i].lhs)
                push!(rm_idxs, i)
                break
            end
        end
    end
    deleteat!(eqs, rm_idxs)

    # draw tethers
    @parameters begin
        m_elevation = elevation
        m_tether_length = tether_length
    end
    @variables begin
        t_x[1:3]
        t_y[1:3]
        t_z[1:3]
        α
        γ[1:se.segments-1]
        local_z[1:se.segments-1]
        local_x[1:se.segments-1]
        h
        r
        total_acc
        kite_distance
    end
    eqs = [
        eqs
        t_z ~ rotate_in_xz([1.0, 0.0, 0.0], m_elevation)
        t_y ~ [0, 1, 0]
        t_x ~ t_y × t_z
        h ~ kite_distance/(2tan(α))
        r ~ kite_distance/(2sin(α))
        pos[:, end] ~ kite_distance * t_z
        pos[:, 1] ~ p1
        total_acc ~ norm(ss.acc)
    ]
    for i in 1:se.segments-1
        eqs = [
            eqs
            γ[i]         ~ -α + 2α * i / (se.segments+1)
            local_z[i]   ~ (kite_distance/2 + r * sin(γ[i]))
            local_x[i]   ~ (r * cos(γ[i]) - h)
            pos[:, i+1] ~ local_z[i] * t_z + local_x[i] * t_x
        ]
    end
    for i in 1:se.segments+1
        eqs = [
            eqs
            ss.vel[:, i] ~ zeros(3)
            ]
        end
    eqs = [eqs; ss.l_spring ~ m_tether_length/se.segments]
    eqs = reduce(vcat, Symbolics.scalarize.(eqs))
    opt_params = [parameters(ss)..., m_elevation, m_tether_length]
    @show opt_params α kite_distance
    @mtkbuild sys = OptimizationSystem(total_acc, [α, kite_distance], opt_params; constraints=eqs) simplify=false
    ps = [ss.rel_compression_stiffness => se.rel_compression_stiffness, ss.end_mass => 1.0, m_elevation => elevation, m_tether_length => tether_length]
    
    prob = OptimizationProblem(sys, [0.0, tether_length], ps; grad = true, hess = true)
    @time sol = solve(prob, Optimization.LBFGS(); maxiters=1000)
    POS0, VEL0 = sol[ss.pos], sol[ss.vel]
    POS0, VEL0, eqs
end

# rotate a 3d vector around the y axis
function rotate_in_xz(vec, angle)
    result = [
        cos(angle) * vec[1] - sin(angle) * vec[3],
        vec[2],
        cos(angle) * vec[3] + sin(angle) * vec[1]
    ]
    result
end

# rotate a 3d vector around the z axis
function rotate_in_yx(vec, angle)
    result = [
        cos(angle) * vec[1] + sin(angle) * vec[2]
        cos(angle) * vec[2] - sin(angle) * vec[1]
        vec[3]
    ]
    result
end

function model(se, p1, fix_p1, fix_p2)
    mass_per_meter = se.rho_tether * π * (se.d_tether/2000.0)^2
    @parameters rel_compression_stiffness = se.rel_compression_stiffness end_mass=1.0
    @variables begin 
        pos(t)[1:3, 1:se.segments+1]
        vel(t)[1:3, 1:se.segments+1]
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
            if ! fix_p2
                push!(eqs, acc[:, i]         ~ se.g_earth + total_force[:, i] / (0.5 * m_tether_particle + end_mass))
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
    eqs = [D(l_spring)   ~ se.v_ro / se.segments,
    c_spring          ~ se.c_spring / l_spring,
    m_tether_particle ~ mass_per_meter * l_spring,
    damping           ~ se.damping  / l_spring]
    eqs2 = vcat(eqs2, reduce(vcat, eqs))  
    eqs2 = reduce(vcat, Symbolics.scalarize.(eqs2))

    @mtkbuild sys = ODESystem(reduce(vcat, Symbolics.scalarize.(eqs2)), t)
    sys, pos, vel, len, c_spr, eqs2
end

function absmax(a, b)
    abs(a) > abs(b) ? a : b
end

function simulate(se, simple_sys, p1, POS0, VEL0)
    global u0map, integ
    u0map = [
        pos => POS0
        vel => VEL0
    ]
    prob = ODEProblem(simple_sys, u0map, tspan; fully_determined=true)
    integ = init(prob, FBDF(autodiff=true); dt, abstol=tol, reltol=tol, saveat=ts)

    toc()
    elapsed_time = @elapsed step!(integ, 1e-3, true)
    elapsed_time = @elapsed step!(integ, se.duration, true)

    integ.sol, elapsed_time
end

function play(se, sol, pos)
    dt = 0.151
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
        plot2d(sol[pos][i], time; segments=se.segments, xlim, ylim, xy)
        if se.save
            ControlPlots.plt.savefig("video/"*"img-"*lpad(j, 4, "0"))
        end
        j += 1
        wait_until(start + 0.5 * time * 1e9)
    end
    if se.save
        println("Run the script ./bin/export_gif to create the gif file!")
    end
    nothing
end

function main(; p1=[0,0,0], elevation=0.0, tether_length=10.0, fix_p1=true, fix_p2=false)
    global sol, pos, vel, len, c_spr, simple_sys, eqs
    se = Settings3()
    set_tether_diameter!(se, se.d_tether) # adapt spring and damping constants to tether diameter
    simple_sys, pos, vel, len, c_spr, eqs = model(se, p1, fix_p1, fix_p2)
    POS0, VEL0, eqs = calc_initial_state(se, simple_sys, eqs, p1, elevation, tether_length)
    @assert false
    sol, elapsed_time = simulate(se, simple_sys, p1, POS0, VEL0)
    if @isdefined __PC
        return sol, pos, vel, simple_sys
    end
    play(se, sol, pos)
    println("Elapsed time: $(elapsed_time) s, speed: $(round(se.duration/elapsed_time)) times real-time")
    println("Number of evaluations per step: ", round(sol.stats.nf/(se.duration/0.02), digits=1))
    sol, pos, vel, simple_sys
end

main(elevation=-deg2rad(30), fix_p2=false);

nothing

