# The equilibrium force of a wind-loaded tether

An analytical formula for the axial force in a tether whose two end points are held a fixed
distance apart while a cross wind blows on it — including the case where that distance is
*shorter* than the unstretched tether, which one would expect to put the tether in
compression.

It reproduces the segmented model of [`TetherComponent.jl`](../src/TetherComponent.jl) to
within **0.6%** over 672 operating points spanning a 30x range of lengths, 4x of diameters
and 3x of wind speeds, with no fitted constant. It is implemented as `analytic_force` in
[`examples/plot_compression.jl`](../examples/plot_compression.jl); the data it is checked
against is produced by [`examples/test_compression.jl`](../examples/test_compression.jl).

## The setup

A tether of `n` segments hangs between two anchors, `L` apart on the vertical axis. Its
unstretched length is `L₀`. The wind is horizontal, perpendicular to the line between the
anchors. Gravity is off, so the wind is the only transverse load, and there is no reel-out.

The quantity of interest is the **mean axial force** `F` over the segments, positive under
tension.

| symbol | meaning | unit |
|---|---|---|
| `L`   | distance between the anchors (`l_tether`) | m |
| `L₀`  | unstretched tether length (`l_tether_unstretched`) | m |
| `r`   | `L₀/L`; `r > 1` means the anchors are closer together than the tether is long | – |
| `n`   | number of segments | – |
| `EA`  | axial stiffness; this is exactly the unit spring constant `se.c_spring` | N |
| `w`   | drag per unit length, `½ ρ c_d d v²` | N/m |
| `s`   | sag, the lateral deflection at mid-span | m |
| `F`   | mean axial force, tension positive | N |

`EA = c_spring` because the model's segment stiffness is `c_spring / l_seg` and its force is
`c_spring · (len - l_seg)/l_seg`, i.e. `c_spring` multiplied by strain. Both `EA` and `w`
scale with the diameter, but differently: `EA ∝ d²` (cross section) and `w ∝ d` (frontal
area). That difference is what produces the three regimes below.

## Why a "compressed" tether is in tension

The result that makes the rest of the derivation possible: **the force never actually
becomes compressive**.

A slack cable cannot carry a transverse load at all. So when the anchors are moved closer
together than the unstretched length, the tether does not sit there compressed — it bows out
until its arc is long enough to be *taut again*, and the wind keeps bowing it until the
tension balances the drag. The arc length therefore ends up slightly longer than `L₀`, every
segment is stretched, and the axial force is positive.

The measured sweep confirms it: over all 672 operating points there is not one with a
non-tensile mean force, and the minimum is 0.068 N — at 2 mm, 10 m/s, `L₀ = 1 m` and 10%
compression, the case that minimises it in every respect at once.

This is why the soft compression branch of the spring
(`rel_compression_stiffness`, 1% of the taut stiffness) never enters the formula: it is
never reached.

## Derivation

Three relations close the system.

**1. Sag, from the transverse force balance.** A cable under a transverse load that is
uniform along its chord hangs in a parabola, with mid-span deflection

```math
s = \frac{w L^2}{8 F}
```

**2. Arc length.** For a parabola of sag `s` over a span `L`, to leading order in `s/L`,

```math
\Delta S_\text{smooth} = S - L = \frac{8 s^2}{3 L}
```

but the model is a chain of straight segments, not a smooth curve. A polyline through the
parabola at `n` equally spaced stations is shorter, and by an amount that can be written
down exactly. With `y = 4s\,t(1-t)` and `t_i = i/n`,

```math
\Delta y_i = \frac{4s}{n^2}\,(n - 1 - 2i), \qquad
\sum_{i=0}^{n-1}(n-1-2i)^2 = \frac{n(n^2-1)}{3}
```

so that

```math
\Delta S = \frac{n}{2L}\sum \Delta y_i^2 = \frac{8 s^2}{3 L}\left(1 - \frac{1}{n^2}\right)
```

The factor `1 - 1/n²` is therefore **derived, not fitted**. For `n = 6` it is `35/36`.
Dropping it leaves a systematic 1.4% error — it is not a detail.

There is no approximation in *which* curve is sampled: a chain of straight segments under
uniform point loads is a funicular polygon, which for a uniform load is exactly a parabola
sampled at the nodes.

**3. The elastic law.** The arc has to be exactly as long as the stretched tether:

```math
\Delta S = L_0\left(1 + \frac{F}{EA}\right) - L
```

**Eliminating `s` and `ΔS`** and dividing by `L` gives, with `r = L₀/L`,

```math
\boxed{\;\frac{r}{EA}F^3 + (r-1)F^2 = \left(1 - \frac{1}{n^2}\right)\frac{w^2 L^2}{24}\;}
```

The linear term vanishes, which is what keeps the closed form short. The left-hand side is
negative at `F = 0` and increases without bound, so there is **exactly one positive root**;
`analytic_force` returns it with Cardano's formula, using the trigonometric branch when the
discriminant says all three roots are real.

## The continuum limit, `n → ∞`

The `1 - 1/n²` is an artefact of the model being a chain of straight segments, not physics.
A real tether is a curve, so the physical formula is the one with the factor gone:

```math
\frac{r}{EA}F^3 + (r-1)F^2 = \frac{w^2 L^2}{24}
```

and the three limits below simplify to

```math
F_\text{ext} = EA\,\frac{1-r}{r}, \qquad
F_{r=1} = \sqrt[3]{\frac{EA\,w^2L^2}{24}}, \qquad
F_\text{comp} = \frac{wL}{\sqrt{24\,(r-1)}}
```

This is the form to use for a physical tether, or for a segmented model with enough
segments. `analytic_force` produces it directly, because `1 - 1/Inf^2` is exactly `1`:

```julia
analytic_force(; v_wind=10.0, d_tether=4.0, l_unstretched=6.6, l_tether=6.0, segments=Inf)
```

**How much the factor is worth.** A model with `n` segments carries less arc length for the
same sag, so it needs more force; using the continuum formula on it *over*-predicts by

```math
\left(1 - \tfrac{1}{n^2}\right)^{-1/2} - 1 \;\approx\; \frac{1}{2n^2}
  \quad\text{(compression)}, \qquad
\left(1 - \tfrac{1}{n^2}\right)^{-1/3} - 1 \;\approx\; \frac{1}{3n^2}
  \quad (r = 1)
```

and not at all in extension, where `n` does not appear. The error dies as `1/n²`:

| `n` | 3 | 4 | 6 | 10 | 20 | 50 |
|---|---|---|---|---|---|---|
| compression | 6.07% | 3.28% | 1.42% | 0.50% | 0.13% | 0.02% |
| `r = 1` | 4.00% | 2.17% | 0.94% | 0.34% | 0.08% | 0.01% |

Against the 6-segment sweep the continuum formula gives a median error of 0.94% and a worst
case of 1.41%, all of it one-sided over-prediction — exactly the `1/(2n²)` bias above. Keep
the factor when comparing with a segmented model; drop it when you want the physics. From
about 20 segments on, the distinction stops mattering.

## The three regimes

Each limit of the same cubic is a regime with its own diameter scaling, and all three are
visible in the measured data:

| regime | which term drops | force (`n` segments) | scales as |
|---|---|---|---|
| **extension**, `r < 1` | drag, against a stiff spring | `F = EA (1-r)/r` | `d²` |
| **zero strain**, `r = 1` | the quadratic term | `F = ((1-1/n²)·EA·w²L²/24)^(1/3)` | `d⁴ᐟ³` |
| **compression**, `r > 1` | the elastic term | `F = wL / √(24(r-1)/(1-1/n²))` | `d¹` |

(For a real tether drop the `1 - 1/n²`, as above; it does not change any exponent.)

Extension is pure Hooke's law: the tether is pulled taut, the drag is negligible against the
spring, and the force follows the cross section. Compression is pure drag: the stiffness has
dropped out entirely and the force is set by how hard the wind pushes on the frontal area.
The crossover at `r = 1` is the geometric mean of the two — `EA ∝ d²` and `w² ∝ d²` under a
cube root give `d⁴ᐟ³`.

The measured exponents, taken by dividing each diameter's force by the 2 mm value at fixed
wind, length and ratio:

| ratio | 4 mm / 2 mm | 6 mm / 2 mm | 8 mm / 2 mm | exponent |
|---|---|---|---|---|
| 0.998 (extension) | 4.00 | 9.00 | 16.00 | `d²` |
| 1.000 (zero strain) | 2.52 | 4.33 | 6.35 | `d⁴ᐟ³` |
| ≥ 1.02 (compression) | 2.00 | 3.00 | 4.00 | `d¹` |

`2^(4/3) = 2.5198`, `3^(4/3) = 4.3267`, `4^(4/3) = 6.3496`.

The compression branch also gives the rest of the scaling, all of it confirmed by the
sweep: `F ∝ v²` (measured ratios of exactly 4.00 and 9.00 going from 10 to 20 and 30 m/s),
`F ∝ L`, and `F ∝ 1/√(r-1)`.

## Why the sign never changes

The question this formula was written to answer. In the compression branch,

```math
F \;\to\; \frac{w L}{\sqrt{24 (r-1)/(1 - 1/n^2)}}
```

which is strictly positive for any `r > 1` and any non-zero `w`. The force reaches zero only
in the limit `w → 0` — that is, only when the cross wind vanishes. There is no combination
of length, diameter, wind speed or compression at which it changes sign; it merely gets
small, as `1/√(r-1)` flattens out and `w` shrinks.

## Accuracy

Against the 672 measured operating points (`L₀` ∈ {1, 3, 10, 30} m, `d` ∈ {2, 4, 6, 8} mm,
`v` ∈ {10, 20, 30} m/s, `r` ∈ 0.99 … 1.10):

| | relative error |
|---|---|
| median | 0.003% |
| 90th percentile | 0.237% |
| worst case | 0.606% |

Two spot checks, both to the digits printed:

- extension, `r = 0.998`, `d = 2 mm`: `EA(1-r)/r = 307.92 N` against a measured 307.9167 N.
- zero strain, `r = 1`, `L₀ = 30 m`, 10 m/s: the closed form gives 42.571 / 107.271 /
  184.193 / 270.307 N for 2 / 4 / 6 / 8 mm, against measured 42.571 / 107.272 / 184.193 /
  270.308 N.

The residual is one-sided (`-0.61% … 0.00%`), so it is a systematic shortfall of the
approximation rather than scatter, and all of the worst cases sit at `r = 1.10`, the deepest
compression, where the sag reaches about 20% of the span.

## Assumptions and limits

The formula is exact only in the limit of small sag and uniform transverse load. What breaks
first, in order:

- **Small sag.** Both the parabolic shape and the `8s²/3L` arc length are leading order in
  `s/L`. At `r = 1.10` the sag is ~20% of the span, which is where the 0.6% lives. Using the
  exact parabola arc length does *not* help — it removes a cancellation with the next point
  and makes things worse (6.3%).
- **Uniform load.** The model's drag uses the wind component perpendicular to each segment,
  which falls off as the end segments tilt away from vertical, so `w` is not quite uniform
  along a deeply bowed tether.
- **One `F` for two roles.** The sag relation wants the chord component of the force, the
  elastic relation wants the mean tension along the arc. They coincide to leading order, and
  empirically the single `F` of the cubic tracks the **mean** segment force to 0.6% — but it
  is *not* the anchor force, which it misses by up to 32% at high sag. Compare against
  `f_mean`, not `f_top`/`f_bot`.
- **No gravity, no reel-out, wind perpendicular to the tether**, and all segments taut. The
  last one holds automatically for `r > 1`, as shown above, but it is an assumption for
  `r < 1` with a very large sag.

## Usage

```julia
include("examples/plot_compression.jl")

analytic_force(; v_wind=10.0, d_tether=4.0, l_unstretched=6.6, l_tether=6.0, segments=6)
```

To reproduce the validation, run the sweep and check the formula against it:

```julia
include("examples/test_compression.jl")   # ~3 minutes, writes the CSV, then plots
```

or, on a CSV that already exists,

```julia
include("examples/plot_compression.jl")   # read, check, plot; no ModelingToolkit needed
```

`check_formula(results)` prints the error table above and `plot_lengths(results)` draws the
formula as a dashed line over every measured curve.
