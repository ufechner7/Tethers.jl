__BENCH__ = true; 
let
    include("Tether_06c.jl")
    set = deepcopy(Settings2())
    set.duration = 10.0
    set.callbacks = false
    simple_sys, pos, vel = model(set)
    sol, elapsed_time = simulate(set, simple_sys)
    println("Tethers_06c, without callbacks")
    println("Elapsed time: $(elapsed_time) s, speed: $(round(set.duration/elapsed_time)) times real-time\n")
end

__BENCH__ = true; 
let
    include("Tether_07.jl")
    set = deepcopy(Settings3())
    set.duration = 10.0
    simple_sys, pos, vel = model(set)
    sol, elapsed_time = simulate(set, simple_sys)
    println("Tethers_07, with tether drag")
    println("Elapsed time: $(elapsed_time) s, speed: $(round(set.duration/elapsed_time)) times real-time")
end