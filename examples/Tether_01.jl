# Example one: Falling mass.
using Timers
tic()
using ModelingToolkit, OrdinaryDiffEq, MakieControlPlots
using ModelingToolkit: t_nounits as t, D_nounits as D
using Tethers: display_if_interactive
toc()

G_EARTH::Vector{Float64} = [0.0, 0.0, -9.81]    # gravitational acceleration     [m/s²]

# defining the model
@variables   pos(t)[1:3] = [0.0, 0.0,  0.0]
@variables   vel(t)[1:3] = [0.0, 0.0, 50.0] 
@variables   acc(t)[1:3]

eqs = vcat(D(pos) ~ vel,
           D(vel) ~ acc,
           acc    ~ G_EARTH)

@named sys = System(eqs, t)
simple_sys = mtkcompile(sys)

# running the simulation
duration = 10.0
dt       = 0.02
tol      = 1e-6
ts       = 0:dt:duration

prob = ODEProblem(simple_sys, nothing, (0.0, duration))
toc()
@time sol = solve(prob, Rodas5(), dt=dt, abstol=tol, reltol=tol, saveat=ts)

# plotting the result
X = sol.t
POS_Z = stack(sol[pos], dims=1)[:,3]
VEL_Z = stack(sol[vel], dims=1)[:,3]

p = plot(X, POS_Z, VEL_Z; xlabel="time [s]", ylabels=["pos_z [m]", "vel_z [m/s]"],
         labels=["pos_z [m]", "vel_z [m/s]"], fig="falling mass", xticks=0:2:duration,
         yticks=20)
display_if_interactive(p)

# saving the result for comparison with the Python implementation
mkpath("output")
open(joinpath("output", "Tether_01_julia.csv"), "w") do io
    println(io, "time,pos_z,vel_z")
    for i in eachindex(X)
        println(io, "$(X[i]),$(POS_Z[i]),$(VEL_Z[i])")
    end
end
