# Quasi-steady tether shape for the MATLAB reference case, compared against the analytic
# catenary that the same case produces without drag.
#
# `import`, not `using`: menu.jl runs every example into the same `Main`, and GLMakie
# exports `plot` just like MakieControlPlots, so a `using GLMakie` here would make `plot`
# ambiguous in every example run afterwards.
import GLMakie
using Tethers: display_if_interactive
using Tethers.QuasiSteady: get_initial_conditions, get_analytic_catenary, Tether, step!

const DATA = joinpath(@__DIR__, "..", "..", "test", "data")

# Read the initial conditions from the .mat file. The fixture already supplies a
# state_vec/kite_pos/wind_vel of its own, so there is nothing left for `init!` to derive -
# wrap the `settings` in a `Tether`, seed its state directly, and solve with one `step!`.
state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings = get_initial_conditions(joinpath(DATA, "input_basic_test.mat"))
te = Tether(settings)
te.wind_vel .= wind_vel
te.state_vec .= state_vec
step!(te, kite_pos, kite_vel; tether_length, wind_vel)

x_qs = collect(sqrt.(te.tether_pos[1,:].^2 + te.tether_pos[2,:].^2))
y_qs = collect(te.tether_pos[3,:])

# Read the catenary curve
x_cat, y_cat = get_analytic_catenary(joinpath(DATA, "input_analytic_catenary.mat"))

fig = GLMakie.Figure()
ax1 = GLMakie.Axis(fig[1, 1]; xlabel="X [m]", ylabel="Y [m]", autolimitaspect=1)
l_cat = GLMakie.lines!(ax1, x_cat, y_cat)
l_qs  = GLMakie.scatterlines!(ax1, x_qs, y_qs)
GLMakie.Legend(fig[1, 2], [l_cat, l_qs], ["Analytic catenary", "Quasi-steady model"])

tether_pos = hcat(te.p0, te.tether_pos, [0; 0; 0])
ax2 = GLMakie.Axis3(fig[2, 1]; title="3D view", xlabel="X [m]", ylabel="Y [m]", zlabel="Z [m]",
                    aspect=:data)
s_tether = GLMakie.scatterlines!(ax2, tether_pos[1,:], tether_pos[2,:], tether_pos[3,:])
s_origin = GLMakie.scatter!(ax2, [0.0], [0.0], [0.0]; markersize=20, marker=:rect, color=:gray)
s_kite   = GLMakie.scatter!(ax2, [te.p0[1]], [te.p0[2]], [te.p0[3]]; markersize=12, marker=:diamond, color=:green)
GLMakie.Legend(fig[2, 2], [s_tether, s_origin, s_kite], ["Tether", "Origin", "Kite"])
display_if_interactive(fig)
