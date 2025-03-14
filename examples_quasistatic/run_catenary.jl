include("../src/Tether_quasistatic.jl")

# Read the initial conditions from a .mat file
state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings = get_initial_conditions("test/data/input_basic_test.mat")
println(state_vec)
state_vec, tether_pos, Ft_ground, Ft_kite, p0 =  simulate_tether(state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings)
println(state_vec)

x_qs = vec(sqrt.(tether_pos[1,:].^2 + tether_pos[2,:].^2))
y_qs = vec(tether_pos[3,:])

# Read the catenary curve
x_cat, y_cat = get_analytic_catenary("test/data/input_analytic_catenary.mat")

plt.plot(x_cat, y_cat)
plt.plot(x_qs, y_qs)
plt.legend(["Analytic catenary", "Quasi static model"])
plt.grid(true)
plt.show()