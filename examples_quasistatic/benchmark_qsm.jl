using BenchmarkTools

include("../src/Tether_quasistatic.jl")

# Read the initial conditions from a .mat file
state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings = get_initial_conditions("test/data/input_basic_test.mat")

@benchmark simulate_tether(state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings)






