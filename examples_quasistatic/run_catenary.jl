include("../src/Tether_quasistatic.jl")

# Get information from the .mat
state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings = get_initial_conditions("test/data/input_basic_test.mat")
Fobj_ref, p0_ref, pj_ref, T0_ref = get_test_output("test/data/basic_test_results.mat")

state_vec, tether_pos, Ft_ground, Ft_kite, p0 =  simulate_tether(state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings)

x_qs = vec(sqrt.(tether_pos[1,:].^2 + tether_pos[2,:].^2))
y_qs = vec(tether_pos[3,:])

# Get catenary curve
x_cat, y_cat = get_analytic_catenary("test/data/input_analytic_catenary.mat")

plt.plot(x_cat, y_cat)
plt.plot(x_qs, y_qs)
plt.legend(["Analytic catenary", "Quasi static model"])
plt.grid(true)
plt.show()