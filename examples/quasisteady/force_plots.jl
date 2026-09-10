# Tether shape and tether force as a function of the kite distance, from the quasi-steady
# model.
#
# `import`, not `using`: menu.jl runs every example into the same `Main`, and GLMakie
# exports `plot` just like MakieControlPlots, so a `using GLMakie` here would make `plot`
# ambiguous in every example run afterwards.
using StaticArrays, LinearAlgebra
import GLMakie
using Tethers: display_if_interactive
using Tethers.QuasiSteady: get_initial_conditions, Tether, step!

function main()
    # Read the initial conditions from a .mat file. The fixture supplies its own settings
    # (including segments) and initial state_vec, so there is nothing left for `init!` to
    # derive - seed a `Tether` directly and solve with `step!`.
    state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings =
        get_initial_conditions(joinpath(@__DIR__, "..", "..", "test", "data", "input_basic_test.mat"))
    te = Tether(settings)
    te.wind_vel .= wind_vel
    te.state_vec .= state_vec

    kite_pos = MVector(5, 100, 300)
    step!(te, kite_pos, kite_vel; tether_length, wind_vel)
    tether_pos = Matrix(te.tether_pos)

    fig1 = GLMakie.Figure()
    ax = GLMakie.Axis3(fig1[1, 1]; title="3D view", xlabel="X [m]", ylabel="Y [m]", zlabel="Z [m]")
    l_tether = GLMakie.scatterlines!(ax, tether_pos[1,:], tether_pos[2,:], tether_pos[3,:])
    s_origin = GLMakie.scatter!(ax, [0.0], [0.0], [0.0]; markersize=20, marker=:rect, color=:gray)
    s_kite   = GLMakie.scatter!(ax, [te.p0[1]], [te.p0[2]], [te.p0[3]]; markersize=12, marker=:diamond, color=:green)
    GLMakie.Legend(fig1[1, 2], [l_tether, s_origin, s_kite], ["Tether", "Origin", "Kite"])
    display_if_interactive(fig1)

    fig2 = GLMakie.Figure()
    ax = GLMakie.Axis(fig2[1, 1]; title="2D view", xlabel="X [m]", ylabel="Z [m]")
    GLMakie.scatterlines!(ax, sqrt.(tether_pos[1,:].^2 + tether_pos[2,:].^2), tether_pos[3,:])
    GLMakie.scatter!(ax, [0.0], [0.0]; markersize=20, marker=:rect, color=:gray)
    GLMakie.scatter!(ax, [sqrt(te.p0[1]^2 + te.p0[2]^2)], [te.p0[3]]; markersize=12, marker=:diamond, color=:green)
    display_if_interactive(fig2)

    x_positions = LinRange(100, 300, 6)

    fig3 = GLMakie.Figure()
    ax1 = GLMakie.Axis(fig3[1, 1]; title="2D view for different kite positions", xlabel="X [m]", ylabel="Z [m]")
    fig4 = GLMakie.Figure()
    ax2 = GLMakie.Axis(fig4[1, 1]; title="Force vs. kite distance", xlabel="Distance [m]", ylabel="Force [N]")

    local l_tether1, s_ground1, s_gnd_force, s_kite_force
    for ii = 1:length(x_positions)
        kite_pos = MVector(5, x_positions[ii], 300)
        step!(te, kite_pos, kite_vel; tether_length, wind_vel)
        tether_pos = Matrix(te.tether_pos)

        l_tether1 = GLMakie.scatterlines!(ax1, sqrt.(tether_pos[1,:].^2 + tether_pos[2,:].^2), tether_pos[3,:]; color=:blue)
        s_ground1 = GLMakie.scatter!(ax1, [0.0], [0.0]; markersize=20, marker=:rect, color=:gray)
        GLMakie.scatter!(ax1, [sqrt(te.p0[1]^2 + te.p0[2]^2)], [te.p0[3]]; markersize=12, marker=:diamond, color=:green)

        distance = sqrt(te.p0[1]^2 + te.p0[2]^2 + te.p0[3]^2)
        s_gnd_force  = GLMakie.scatter!(ax2, [distance], [te.force_gnd]; color=(:blue, 0.5), markersize=10)
        s_kite_force = GLMakie.scatter!(ax2, [distance], [norm(te.force_kite)]; color=(:orange, 0.5), markersize=10)
    end
    GLMakie.Legend(fig3[1, 2], [l_tether1, s_ground1], ["Tether", "Ground station"])
    GLMakie.Legend(fig4[1, 2], [s_gnd_force, s_kite_force], ["Tension at ground station", "Tension at kite"])
    display_if_interactive(fig3)
    display_if_interactive(fig4)
    nothing
end
