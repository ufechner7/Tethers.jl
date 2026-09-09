# Quasi-static tether shape for the MATLAB reference case, compared against the analytic
# catenary that the same case produces without drag.
#
# `import`, not `using`: menu.jl runs every example into the same `Main`, and GLMakie
# exports `plot` just like MakieControlPlots, so a `using GLMakie` here would make `plot`
# ambiguous in every example run afterwards.
import GLMakie
using Tethers: display_if_interactive
include("../../src/Tether_quasistatic.jl")

const DATA = joinpath(@__DIR__, "..", "..", "test", "data")

# Read the initial conditions from the .mat file
state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings = get_initial_conditions(joinpath(DATA, "input_basic_test.mat"))
state_vec, tether_pos, Ft_ground, Ft_kite, p0 =  simulate_tether(state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings)


x_qs = vec(sqrt.(tether_pos[1,:].^2 + tether_pos[2,:].^2))
y_qs = vec(tether_pos[3,:])

# Read the catenary curve
x_cat, y_cat = get_analytic_catenary(joinpath(DATA, "input_analytic_catenary.mat"))

fig1 = GLMakie.Figure()
ax1 = GLMakie.Axis(fig1[1, 1]; xlabel="X [m]", ylabel="Y [m]", autolimitaspect=1)
l_cat = GLMakie.lines!(ax1, x_cat, y_cat)
l_qs  = GLMakie.scatterlines!(ax1, x_qs, y_qs)
GLMakie.Legend(fig1[1, 2], [l_cat, l_qs], ["Analytic catenary", "Quasi static model"])
display_if_interactive(fig1)


tether_pos = hcat(p0, tether_pos, [0; 0; 0])
fig2 = GLMakie.Figure()
ax2 = GLMakie.Axis3(fig2[1, 1]; title="3D view", xlabel="X [m]", ylabel="Y [m]", zlabel="Z [m]",
                    aspect=:data, limits=((0, 100), (0, 100), (0, 800)))
s_tether = GLMakie.scatter!(ax2, tether_pos[1,:], tether_pos[2,:], tether_pos[3,:])
s_origin = GLMakie.scatter!(ax2, [0.0], [0.0], [0.0]; markersize=20, marker=:rect, color=:gray)
s_kite   = GLMakie.scatter!(ax2, [p0[1]], [p0[2]], [p0[3]]; markersize=12, marker=:diamond, color=:green)
GLMakie.Legend(fig2[1, 2], [s_tether, s_origin, s_kite], ["Tether", "Origin", "Kite"])
display_if_interactive(fig2)
