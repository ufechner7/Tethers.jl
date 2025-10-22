using BenchmarkTools

include("../src/Tether_qsm_dual.jl")

# Read the initial conditions from a .mat file
state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings = get_initial_conditions("test/data/input_basic_test.mat")
simulate_tether(state_vec, kite_pos, kite_vel, MMatrix{3, 15}(wind_vel), tether_length, settings; prn=true)

@benchmark simulate_tether(state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings)

#= On Ryzen 7950X 
Iterations: 36
BenchmarkTools.Trial: 10000 samples with 1 evaluation per sample.
 Range (min … max):  103.920 μs …   2.508 ms  ┊ GC (min … max):  0.00% … 91.64%
 Time  (median):     124.560 μs               ┊ GC (median):     0.00%
 Time  (mean ± σ):   149.064 μs ± 169.569 μs  ┊ GC (mean ± σ):  12.93% ± 10.44%

  █▂                                                             
  ██▂▄▂▂▂▂▂▁▂▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▂▂▂▂▂▂▂ ▂
  104 μs           Histogram: frequency by time         1.49 ms <

 Memory estimate: 328.67 KiB, allocs estimate: 2292.
 =#






