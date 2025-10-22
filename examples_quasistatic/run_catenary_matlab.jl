using ControlPlots
include("../src/Tether_quasistatic.jl")

# Read the initial conditions from the .mat file
state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings = get_initial_conditions("test/data/input_basic_test.mat")
state_vec, tether_pos, Ft_ground, Ft_kite, p0 =  simulate_tether(state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings)


x_qs = vec(sqrt.(tether_pos[1,:].^2 + tether_pos[2,:].^2))
y_qs = vec(tether_pos[3,:])

# Read the catenary curve
x_cat, y_cat = get_analytic_catenary("test/data/input_analytic_catenary.mat")

fig, ax = plt.subplots(1,1)
ax.plot(x_cat, y_cat)
ax.plot(x_qs, y_qs, marker = "o")
ax.set_aspect("equal")
ax.legend(["Analytic catenary", "Quasi static model"])
plt.xlabel("X [m]")
plt.ylabel("Y [m]")
fig.show()


tether_pos = hcat(p0, tether_pos, [0; 0; 0])
plt.figure("3D view").add_subplot(projection="3d").set_aspect("equal")
plt.scatter3D(tether_pos[1,:], tether_pos[2,:], tether_pos[3,:], marker = "o")
plt.scatter3D(0, 0, 0,  s = 200, marker = "s", c = "C7")
plt.scatter3D(p0[1], p0[2], p0[3], s = 50, marker = "D", c = "g")
plt.xlabel("X [m]")
plt.ylabel("Y [m]")
plt.xlim(0, 100)
plt.ylim(0, 100)
plt.zlim(0, 800)
plt.legend(["Tether", "Origin", "Kite"])
plt.show()

