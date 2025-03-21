include("../src/Tether_quasistatic.jl")

function main()
    avg_el = deg2rad(70)
    cone_ang = deg2rad(10)
    traj_dist = 500
    gamma = LinRange(0, 2*pi, 20)

    traj_x = traj_dist*sin(cone_ang).*cos.(gamma)
    traj_y = traj_dist*sin(cone_ang).*sin.(gamma)
    traj_z = traj_dist*cos(cone_ang) .+ 0.0.*gamma
    traj = [traj_x'; traj_y'; traj_z']

    rot_mat = [1 0 0; 0 cos(avg_el) -sin(avg_el); 0 sin(avg_el) cos(avg_el)] 

    for ii = 1:size(traj)[2]
    traj[:,ii] .= rot_mat*traj[:,ii]
    end


    fig = plt.figure("3D view")
    ax = fig.add_subplot(111, projection="3d")

    # Scatter plot of the trajectory
    ax.scatter3D(traj[1, :], traj[2, :], traj[3, :])

    # Plot a line from origin to the last point
    #ax.plot3D([0, traj[1, end]], [0, traj[2, end]], [0, traj[3, end]])
    ax.plot3D([0, 0], [0, -traj_dist*sin(avg_el)], [0, traj_dist*cos(avg_el)], ":")
    plt.xlabel("X [m]")
    plt.ylabel("Y [m]")
    plt.zlabel("Z [m]")
    plt.xlim((-250, 250))
    plt.zlim((-50, 450))
    plt.title("Circular trajectory")
    plt.legend(["Kite trajectory", "Average elevation axis"])
    # Show the plot
    plt.show()
    #traj = rot_mat*traj

    # Initial position gamma = 0
    kite_pos = MVector{3}(traj[:, 1])
    # Initialise model
    tether_length = 400
    segments = 10
    state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings = init_quasistatic(kite_pos, tether_length, segments = segments)
    state_vec, tether_pos, Ft_ground, Ft_kite, p0 =  simulate_tether(state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings)

    tether_pos = hcat(p0, tether_pos, [0; 0; 0])
    plt.figure("3D view").add_subplot(projection="3d").set_aspect("equal")
    plt.plot3D(tether_pos[1,:], tether_pos[2,:], tether_pos[3,:], marker = "o")
    plt.scatter3D(0, 0, 0,  s = 200, marker = "s", c = "C7")
    plt.scatter3D(p0[1], p0[2], p0[3], s = 50, marker = "D", c = "g")
    plt.xlabel("X [m]")
    plt.ylabel("Y [m]")
    plt.legend(["Tether", "Origin", "Kite"])
    plt.show()

    all_tether_pos = zeros(length(gamma), 3, segments)
    Ft_ground = zeros(length(gamma), 1)
    
    for ii = 1:length(gamma)
        state_vec, all_tether_pos[ii, :, :], Ft_ground[ii, :], Ft_kite, p0 = simulate_tether(state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings)
    end

    

end