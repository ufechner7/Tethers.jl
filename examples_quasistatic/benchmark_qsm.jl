using BenchmarkTools

include("../src/Tether_quasistatic.jl")

# Get information from the .mat
state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings = getInitCond("test/data/input_basic_test.mat")
Fobj_ref, p0_ref, pj_ref, T0_ref = getOutputObjFun("test/data/basic_test_results.mat")

@benchmark tetherQuasiStatic(state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings)

@time tetherQuasiStatic(state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings) 





