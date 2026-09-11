# Quasi-steady tether shape for a kite at a fixed position, plotted in 3D.
#
# `import`, not `using`: menu.jl runs every example into the same `Main`, and GLMakie
# exports `plot` just like MakieControlPlots, so a `using GLMakie` here would make `plot`
# ambiguous in every example run afterwards.
using StaticArrays, LinearAlgebra
import GLMakie
using Tethers: display_if_interactive
using Tethers.QuasiSteady: StaticSettings, Tether, init!

# Set initial conditions: StaticSettings takes elevation/azimuth/l_tether, KiteUtils-style,
# so the fixed kite position is converted to that form.
kite_pos = MVector{3}([100.0, 100, 800])
β0 = asin(kite_pos[3] / norm(kite_pos))
φ0 = atan(kite_pos[2], kite_pos[1])
se = StaticSettings(segments = 22, elevation = rad2deg(β0), azimuth = rad2deg(φ0),
                    l_tether = 1.05 * norm(kite_pos))
te = Tether(se)
init!(te)

x_qs = vec(sqrt.(te.tether_pos[1,:].^2 + te.tether_pos[2,:].^2))
y_qs = vec(te.tether_pos[3,:])

tether_pos = hcat(te.p0, te.tether_pos, [0; 0; 0])
x_min, x_max = extrema(tether_pos[1,:])
y_min, y_max = extrema(tether_pos[2,:])
z_min, z_max = extrema(tether_pos[3,:])
fig = GLMakie.Figure()
ax = GLMakie.Axis3(fig[1, 1]; title="3D view", xlabel="X [m]", ylabel="Y [m]", zlabel="Z [m]",
                   aspect=:data, limits=((min(0, x_min), max(100, x_max)),
                                          (min(0, y_min), max(100, y_max)),
                                          (min(0, z_min), max(800, z_max))))
l_tether = GLMakie.scatterlines!(ax, tether_pos[1,:], tether_pos[2,:], tether_pos[3,:])
s_origin = GLMakie.scatter!(ax, [0.0], [0.0], [0.0]; markersize=20, marker=:rect, color=:gray)
s_kite   = GLMakie.scatter!(ax, [te.p0[1]], [te.p0[2]], [te.p0[3]]; markersize=12, marker=:diamond, color=:green)
GLMakie.Legend(fig[1, 2], [l_tether, s_origin, s_kite], ["Tether", "Origin", "Kite"])
display_if_interactive(fig)
