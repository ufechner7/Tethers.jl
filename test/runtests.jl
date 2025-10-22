using Test, LinearAlgebra

__BENCH__ = true; 

# @testset "Tether_06c" begin   
#     # without callbacks
#     include("../src/Tether_06c.jl")
#     set = deepcopy(Settings2())
#     set.duration = 10.0
#     set.callbacks = false
#     simple_sys, pos, vel = model(set)
#     sol, elapsed_time = simulate(set, simple_sys)
#     @test elapsed_time < 1.0
#     l_tether_theoretical = set.l0 + set.v_ro * set.duration
#     @test l_tether(sol, pos) ≈ l_tether_theoretical rtol=2e-3
#     events = Int64(round(length(sol.t)- set.duration/set.dt)-1)
#     @test events == 0
    
#     # with callbacks
#     set.callbacks = true
#     simple_sys, pos, vel = model(set)
#     sol, elapsed_time = simulate(set, simple_sys)
#     @test elapsed_time < 1.0
#     l_tether_theoretical = set.l0 + set.v_ro * set.duration
#     @test l_tether(sol, pos) ≈ l_tether_theoretical rtol=2e-3
#     events = Int64(round(length(sol.t)- set.duration/set.dt)-1)
#     @test events >= 4 # 8 events with Rodas5, 4 events with KenCarp4
# end
include("test_qsm.jl")
nothing

