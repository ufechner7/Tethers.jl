using LaTeXStrings
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
    
    # Initial position gamma = 0
    kite_pos = MVector{3}(traj[:, 1])
    # Initialise model
    tether_length = 100
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

    all_tether_pos = zeros(length(gamma), 3, segments + 1)
    all_Ft_kite = zeros(3, length(gamma))
    all_Ft_ground = zeros(length(gamma))
    
    for ii = 1:length(gamma)
        kite_pos .= MVector{3}(traj[:, ii])
        state_vec, tether_pos, Ft_ground, Ft_kite, p0 = simulate_tether(state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings)
        state_vec .= state_vec
        tether_pos = hcat(p0, tether_pos)
        all_tether_pos[ii, :, :] .= tether_pos   
        all_Ft_kite[:, ii] .= Ft_kite
        all_Ft_ground .= Ft_ground
    end
    
    plt.figure("2D plot")
    plt.plot(gamma, all_Ft_kite[1, :]./1000)
    plt.plot(gamma, all_Ft_kite[2, :]./1000)
    plt.plot(gamma, all_Ft_kite[3, :]./1000)
    plt.legend([L"F_x", L"F_y", L"F_z"])    
    plt.xlabel(L"\gamma [rad]")
    plt.ylabel("Force [kN]")
    plt.title("Tether force components at kite during a circular trajectory")
    plt.show()
    

    plt.figure("3D view").add_subplot(projection="3d").set_aspect("equal")
    plt.scatter3D(0, 0, 0,  s = 200, marker = "s", c = "C7")
    plt.scatter3D(traj[1, :], traj[2, :], traj[3, :])
    for ii = 1:length(gamma)
        plt.plot3D(all_tether_pos[ii,1,:], all_tether_pos[ii,2,:], all_tether_pos[ii,3,:], marker = "x", color = "C1", linestyle = ":")
    end
    plt.xlabel("X [m]")
    plt.ylabel("Y [m]")
    plt.legend(["Origin", "Kite trajectory", "Tethers"])
    plt.show()   
end