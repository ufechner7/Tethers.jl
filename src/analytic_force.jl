"""
    analytic_force(se; v_wind_perp, d_segment, l_unstretched, l_segment, segments=Inf)

Step two of PlanCompression.md: the analytical prediction of the mean axial force [N],
tension positive, without solving the model.

`se` is a [`TetherSettings`](@ref) (or anything with the same `rho`, `cd_tether`, `d_tether`
and `c_spring` fields); `rho` and `cd_tether` are taken from it directly, and the axial
stiffness is scaled from `se.c_spring` to the given `d_segment` the same way
[`set_diameter!`](@ref) does, `EA = se.c_spring * (d_segment/se.d_tether)^2`, so `se` need
not already be set to that diameter.

A tether under a transverse load that is uniform along the chord hangs in a parabola, and a
chain of straight segments under uniform point loads is exactly that parabola sampled at
its nodes. Three relations close the system, with `L` the distance between the anchors,
`L0` the unstretched length, `EA` the axial stiffness and `w` the drag per unit length:

- sag from the force balance,        `s  = w L² / (8F)`
- arc length of the sampled parabola `ΔS = (1 - 1/n²) · 8s²/(3L)`
- and the elastic law,               `ΔS = L0 (1 + F/EA) - L`

Eliminating `s` and `ΔS` leaves a cubic in `F` with exactly one positive root, solved here
in closed form:

    (r/EA) F³ + (r-1) F² = (1 - 1/n²) w² L² / 24,   r = L0/L

The `1 - 1/n²` is the only discretisation term: an `n`-segment polyline through a parabola
is that much shorter than the smooth curve, so it needs that much more sag — and hence less
force — to take up the same slack. It is derived, not fitted.

The default `segments=Inf` gives the continuum limit, i.e. the formula for a real tether
rather than for a chain of `n` segments — `1 - 1/Inf^2` is exactly `1`, so the term simply
drops out. Pass the actual segment count to match a segmented model instead: a 6-segment
model over-predicts the continuum limit by about `1/(2n²)`, 1.4%; by 20 segments it is
0.13%.

Reproduces the measured force to 0.6% over the whole sweep; see `check_formula` in
`examples/plot_compression.jl`. The derivation is in `docs/segment_force.md`.
"""
function analytic_force(se; v_wind_perp, d_segment, l_unstretched, l_segment, segments=Inf)
    EA = se.c_spring * (d_segment/se.d_tether)^2               # axial stiffness    [N]
    w  = 0.5 * se.rho * se.cd_tether * (d_segment/1000) * v_wind_perp^2   # drag per meter [N/m]
    r  = l_unstretched / l_segment
    cn = 1 - 1/segments^2
    # F³ + a₂F² + a₀ = 0; a₁ vanishes, which is what makes the closed form short
    a2 = EA * (r - 1) / r
    a0 = -cn * EA * w^2 * l_segment^2 / (24r)
    p  = -a2^2/3
    q  = 2a2^3/27 + a0
    D  = q^2/4 + p^3/27
    if D > 0            # one real root
        return cbrt(-q/2 + sqrt(D)) + cbrt(-q/2 - sqrt(D)) - a2/3
    end                 # three real roots, the physical one is the positive one
    m  = 2sqrt(-p/3)
    th = acos(clamp(3q/(p*m), -1, 1))/3
    maximum(m*cos(th - 2π*k/3) - a2/3 for k in 0:2)
end
