using StaticArrays, LinearAlgebra
include("../src/Tether_quasistatic.jl")

function main()
    # Read the initial conditions from a .mat file
    state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings = get_initial_conditions("test/data/input_basic_test.mat")

    kite_pos = MVector(5, 100, 300)
    state_vec, tether_pos, Ft_ground, Ft_kite, p0 =  simulate_tether(state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings)

    plt.figure("3D view").add_subplot(projection="3d")
    plt.plot3D(tether_pos[1,:], tether_pos[2,:], tether_pos[3,:], marker = "o")
    plt.scatter3D(0, 0, 0,  s = 200, marker = "s", c = "C7")
    plt.scatter3D(p0[1], p0[2], p0[3], s = 50, marker = "D", c = "g")
    plt.legend(["Tether", "Origin", "Kite"])

    plt.figure("2D view").add_subplot()
    plt.plot(sqrt.(tether_pos[1,:].^2 + tether_pos[2,:].^2), (tether_pos[3,:]), marker = "o")
    plt.scatter(0, 0, s = 200, marker = "s", c = "C7")
    plt.scatter(sqrt(p0[1]^2 +  p0[2]^2), p0[3], s = 50, marker = "D", c = "g")
    
    x_positions = LinRange(100, 300, 6)

    ax1 = plt.figure("2D view for different kite positions").add_subplot()
    ax2 = plt.figure("Force vs. kite distance").add_subplot()

    for ii = 1:length(x_positions)
        kite_pos = MVector(5, x_positions[ii], 300)
        state_vec, tether_pos, Ft_ground, Ft_kite, p0 = simulate_tether(state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings)
        
        ax1.plot(sqrt.(tether_pos[1,:].^2 + tether_pos[2,:].^2), (tether_pos[3,:]), marker = "o", c = "C0")
        ax1.scatter(0, 0, s = 200, marker = "s", c = "C7")
        ax1.scatter(sqrt(p0[1]^2 +  p0[2]^2), p0[3], s = 50, marker = "D", c = "g") 

        distance = sqrt(p0[1]^2 + p0[2]^2 + p0[3]^2)
        ax2.scatter(distance, Ft_ground, c = "C0", s=20,  alpha=0.5)
        ax2.scatter(distance, norm(Ft_kite), c = "C1", s=20, alpha=0.5)        
    end
    ax1.legend(["Tether", "Ground station"])
    ax2.legend(["Tension at ground station", "Tension at kite"])
    plt.show()
end

