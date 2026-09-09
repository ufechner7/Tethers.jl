using Pkg
# When running this script directly (not via `examples/menu3.jl`), activate the examples
# environment (one level up) unless it is already active.
if dirname(Pkg.project().path) != normpath(joinpath(@__DIR__, ".."))
    Pkg.activate(joinpath(@__DIR__, ".."))
end
using BenchmarkTools

include("../../src/Tether_qsm_dual.jl")

const segments = 15

# Initial conditions, hardcoded. These are the values that `test/data/input_basic_test.mat`
# used to provide, so the timings below stay comparable; the .mat file itself is only a
# reference fixture for `test/test_qsm.jl` and has no business in a benchmark.
#
# state_vec is (elevation [rad], wind-frame azimuth [rad], ground tension [N]); the two
# angles are the MATLAB-convention pair (0.321750554, -0.306277369) already converted with
# `matlab_to_wind`.
state_vec     = MVector{3}(1.1302856641844843, -0.7853981636973135, 1.60941384e5)
kite_pos      = MVector{3}(100.0, 100.0, 300.0)
kite_vel      = MVector{3}(0.0, 0.0, 0.0)
wind_vel      = zeros(3, segments)
tether_length = 431.66247903554

settings = Settings(1.225,                                # rho          [kg/m³]
                    MVector{3}(0.0, 0.0, -9.8066),        # g_earth      [m/s²]
                    1.2,                                  # cd_tether
                    29.71973504179974,                    # d_tether     [mm]
                    0.6729014779417218,                   # rho_tether
                    8.047069220746364e7)                  # c_spring = E*A, E = 116 GPa, A = 693.7 mm²

simulate_tether(state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings; prn=true)

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






