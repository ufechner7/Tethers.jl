using BenchmarkTools

const segments = 15

include("../src/Tether_quasistatic.jl")

# Read the initial conditions from a .mat file
state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings = get_initial_conditions("test/data/input_basic_test.mat")
simulate_tether(state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings; prn=true)

@benchmark simulate_tether(state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings)

#= On Ryzen 7950X 
julia> include("examples_quasistatic/benchmark_qsm.jl")
Iterations: 36
BenchmarkTools.Trial: 10000 samples with 1 evaluation per sample.
 Range (min … max):   93.927 μs …   8.020 ms  ┊ GC (min … max): 0.00% … 98.07%
 Time  (median):     102.166 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   106.874 μs ± 117.790 μs  ┊ GC (mean ± σ):  2.18% ±  2.32%

     █▇▃▁                                                        
  ▂▃█████▇▇▆▇▆▆▆▆▆▅▅▅▅▅▅▅▅▅▅▅▅▄▅▄▄▄▄▄▄▃▃▃▃▃▃▃▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▁▂ ▄
  93.9 μs          Histogram: frequency by time          133 μs <

 Memory estimate: 46.53 KiB, allocs estimate: 1118. 
 =#






