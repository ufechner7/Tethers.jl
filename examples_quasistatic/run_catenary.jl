include("../src/Tether_quasistatic.jl")

# Read the initial conditions from the .mat file
# state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings
#_, kite_pos, _, _, tether_length, settings = get_initial_conditions("test/data/input_basic_test.mat")
#state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings = get_initial_conditions("test/data/input_basic_test.mat")
# Set initial conditions
kite_pos = MVector{3}([100, 100, 500])
tether_length = 500;
state_vec = MVector{3}([2.827132149347252, -0.8, 64865.86737989248])
state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings = init_quasistatic(kite_pos, tether_length, segments = 20)

state_vec, tether_pos, Ft_ground, Ft_kite, p0 =  simulate_tether(state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings)
println(state_vec)

x_qs = vec(sqrt.(tether_pos[1,:].^2 + tether_pos[2,:].^2))
y_qs = vec(tether_pos[3,:])

# Read the catenary curve
x_cat, y_cat = get_analytic_catenary("test/data/input_analytic_catenary.mat")

#= plt.plot(x_cat, y_cat)
plt.plot(x_qs, y_qs, marker = "o")
plt.legend(["Analytic catenary", "Quasi static model"])
plt.xlabel("X [m]")
plt.ylabel("Y [m]")
plt.show() =#


tether_pos = hcat(p0, tether_pos, [0; 0; 0])
plt.figure("3D view").add_subplot(projection="3d").set_aspect("equal")
plt.plot3D(tether_pos[1,:], tether_pos[2,:], tether_pos[3,:], marker = "o")
plt.scatter3D(0, 0, 0,  s = 200, marker = "s", c = "C7")
plt.scatter3D(p0[1], p0[2], p0[3], s = 50, marker = "D", c = "g")
plt.xlabel("X [m]")
plt.ylabel("Y [m]")
plt.legend(["Tether", "Origin", "Kite"])
plt.show()

