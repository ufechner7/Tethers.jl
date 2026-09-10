# Tether shape and tether force at the kite while the kite flies a circular trajectory.
#
# `import`, not `using`: menu.jl runs every example into the same `Main`, and GLMakie
# exports `plot` just like MakieControlPlots, so a `using GLMakie` here would make `plot`
# ambiguous in every example run afterwards.
using LaTeXStrings, StaticArrays, LinearAlgebra
import GLMakie
using Tethers: display_if_interactive
using Tethers.QuasiSteady: StaticSettings, Tether, init!, step!

# `display(fig)` re-uses the one GLMakie window, so each figure would replace the previous
# one as soon as it is shown; a fresh `Screen` gives every figure a window of its own, and
# `title` names it, because every window is called "Makie" otherwise. The `Screen` is created
# inside the closure so that nothing opens a window on CI.
show_fig(fig, title) = display_if_interactive(() -> display(GLMakie.Screen(; title), fig))

function main()
    avg_el = deg2rad(70)
    cone_ang = deg2rad(10)
    traj_dist = 500
    gamma = LinRange(0, 2*pi, 20)
    gamma_dot = 0.05

    traj_x = traj_dist*sin(cone_ang).*cos.(gamma)
    traj_y = traj_dist*sin(cone_ang).*sin.(gamma)
    traj_z = traj_dist*cos(cone_ang) .+ 0.0.*gamma
    traj   = [traj_x'; traj_y'; traj_z']

    vel_x = traj_dist*sin(cone_ang).*-sin.(gamma)*gamma_dot
    vel_y = traj_dist*sin(cone_ang).*cos.(gamma)*gamma_dot
    vel_z = 0.0.*gamma*gamma_dot
    vel   = [vel_x'; vel_y'; vel_z']

    acc_x = traj_dist*sin(cone_ang).*-cos.(gamma)*gamma_dot # + terms with gamma_ddot = 0
    acc_y = traj_dist*sin(cone_ang).*-sin.(gamma)*gamma_dot # + terms with gamma_ddot = 0
    acc_z = 0.0.*gamma*gamma_dot^2 # + terms with gamma_ddot = 0
    acc   = [acc_x'; acc_y'; acc_z']




    rot_mat = [1 0 0; 0 cos(avg_el) -sin(avg_el); 0 sin(avg_el) cos(avg_el)] 

    for ii = 1:size(traj)[2]
        traj[:,ii] .= rot_mat*traj[:,ii]
        vel[:,ii] .= rot_mat*vel[:,ii]
        acc[:,ii] .= rot_mat*acc[:,ii]
    end
    
    # Initial position gamma = 0
    kite_pos = MVector{3}(traj[:, 1])
    segments = 20

    # Initialize model: StaticSettings carries the initial condition as elevation/azimuth/
    # l_tether, KiteUtils-style, so the trajectory's first point is converted to that form.
    β0 = asin(kite_pos[3] / norm(kite_pos))
    φ0 = atan(kite_pos[2], kite_pos[1])
    se = StaticSettings(segments = segments, elevation = rad2deg(β0), azimuth = rad2deg(φ0),
                        l_tether = 1.05 * norm(kite_pos))
    te = Tether(se)
    init!(te)   # solves the catenary and one step!, leaving te in a consistent, solved state
    tether_pos = hcat(te.p0, te.tether_pos, [0; 0; 0])

    fig1 = GLMakie.Figure()
    ax = GLMakie.Axis3(fig1[1, 1]; title="3D view", xlabel="X [m]", ylabel="Y [m]", zlabel="Z [m]", aspect=:data)
    l_tether = GLMakie.scatterlines!(ax, tether_pos[1,:], tether_pos[2,:], tether_pos[3,:])
    s_origin = GLMakie.scatter!(ax, [0.0], [0.0], [0.0]; markersize=20, marker=:rect, color=:gray)
    s_kite   = GLMakie.scatter!(ax, [te.p0[1]], [te.p0[2]], [te.p0[3]]; markersize=12, marker=:diamond, color=:green)
    GLMakie.Legend(fig1[1, 2], [l_tether, s_origin, s_kite], ["Tether", "Origin", "Kite"])
    show_fig(fig1, "Initial tether shape")

    # `te.tether_pos` holds the `segments - 1` inner nodes, plus the kite point and the
    # ground station added below, which makes `segments + 1` points per tether
    all_tether_pos = zeros(length(gamma), 3, segments + 1)
    all_Ft_kite = zeros(3, length(gamma))
    all_Ft_ground = zeros(length(gamma))

    elapsed_time = @elapsed for ii = 1:length(gamma)
        kite_pos = MVector{3}(traj[:, ii])
        kite_vel = MVector{3}(vel[:, ii])
        step!(te, kite_pos, kite_vel)   # tether_length defaults to (1 + se.slack) * norm(kite_pos)
        tether_pos = hcat(te.p0, te.tether_pos, [0.0; 0.0; 0.0])
        all_tether_pos[ii, :, :] .= tether_pos
        all_Ft_kite[:, ii] .= te.force_kite
        all_Ft_ground[ii] = te.force_gnd
    end
    duration = 2π / gamma_dot   # simulated time spanned by one full revolution
    println("Elapsed time: $(elapsed_time) s, speed: $(round(duration/elapsed_time)) times real-time")

    fig2 = GLMakie.Figure()
    ax = GLMakie.Axis(fig2[1, 1]; title="Tether force components at kite during a circular trajectory",
                      xlabel=L"\gamma [rad]", ylabel="Force [kN]")
    lx = GLMakie.lines!(ax, gamma, all_Ft_kite[1, :]./1000)
    ly = GLMakie.lines!(ax, gamma, all_Ft_kite[2, :]./1000)
    lz = GLMakie.lines!(ax, gamma, all_Ft_kite[3, :]./1000)
    GLMakie.Legend(fig2[1, 2], [lx, ly, lz], [L"F_x", L"F_y", L"F_z"])
    show_fig(fig2, "Tether force at the kite")
    

    fig3 = GLMakie.Figure()
    ax = GLMakie.Axis3(fig3[1, 1]; title="3D view", xlabel="X [m]", ylabel="Y [m]", zlabel="Z [m]", aspect=:data)
    s_origin = GLMakie.scatter!(ax, [0.0], [0.0], [0.0]; markersize=20, marker=:rect, color=:gray)
    s_traj = GLMakie.scatter!(ax, traj[1, :], traj[2, :], traj[3, :])
    local l_tethers
    for ii = 1:length(gamma)
        l_tethers = GLMakie.scatterlines!(ax, all_tether_pos[ii,1,:], all_tether_pos[ii,2,:], all_tether_pos[ii,3,:];
                                          marker=:xcross, color=:orange, linestyle=:dot)
    end
    GLMakie.Legend(fig3[1, 2], [s_origin, s_traj, l_tethers], ["Origin", "Kite trajectory", "Tethers"])
    show_fig(fig3, "Tether shapes along the trajectory")
    nothing
end