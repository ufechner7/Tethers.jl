__BENCH__ = true; include("Tether_07.jl")

set = deepcopy(Settings())
set.duration = 10.0
simple_sys, pos, vel = model(set)
sol, elapsed_time = simulate(set, simple_sys) # warm-up
sol, elapsed_time = simulate(set, simple_sys)

println("Elapsed time: $(elapsed_time) s, speed: $(round(set.duration/elapsed_time)) times real-time")