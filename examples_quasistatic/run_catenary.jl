include("../src/Tether_quasistatic.jl")

# Read the initial conditions from the .mat file
# state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings
#_, kite_pos, _, _, tether_length, settings = get_initial_conditions("test/data/input_basic_test.mat")
state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings = get_initial_conditions("test/data/input_basic_test.mat")
# Set initial conditions
#state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings = init_quasistatic(kite_pos, tether_length, settings=settings, wind_vel=wind_vel)
println(state_vec)
kite_pos = MVector{3}([500, 500, 200])
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

theta = atan(tether_pos[3, end-1]/sqrt(tether_pos[1, end-1]^2 + tether_pos[2, end-1]^2))
phi = atan(tether_pos[2, end-1]/tether_pos[1, end-1])
println("phi - elevation:")
println(phi)
println(rad2deg(phi))

println("Phi:")
println(phi)
println(rad2deg(phi))

plt.figure("3D view").add_subplot(projection="3d")
plt.plot3D(tether_pos[1,:], tether_pos[2,:], tether_pos[3,:], marker = "o")
plt.scatter3D(0, 0, 0,  s = 200, marker = "s", c = "C7")
plt.scatter3D(p0[1], p0[2], p0[3], s = 50, marker = "D", c = "g")
plt.legend(["Tether", "Origin", "Kite"])
plt.show()

