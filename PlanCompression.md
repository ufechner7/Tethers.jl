# Find an analytical formula for the compression force
The compression force will never change its sign as long as the tangential wind is significant.
Try to derive an analytical formula.

## Step one: Investigate
1. create a script examples/test_compression.jl . It shall create a vertical tether with 6 segments and 1m segment length. The wind speed shall be 10 m/s. Turn off gravity. Now measure the equilibrium force for 1% extension to 10% compression and plot it.

2. Can you do the same, but now also vary the tether length, using 1, 3, 10 and 30m as steps
for the unstretched length l_tether_unstretched? And create a result file with the columns
v_wind, l_tether_unstretched, l_tether, f_top, f_bot, f_mean and f_seg_1 … f_seg_6, where
l_tether is the distance between the first and the last point, and vary the ratio
l_tether_unstretched / l_tether from 0.99 (1% extension) to 1.10 (10% compression)?
Furthermore, the plot shall use a logarithmic axis for the force.

3. can you extend the results from step 2. by running the same tests at 20 m/s wind and 30 m/s wind

4. can you extend the results from step 3 by running the same tests with 2,4,6 and 8 mm tether diameter

### Status

The work is split over two scripts, both reachable from `examples/menu.jl`:

- [examples/test_compression.jl](examples/test_compression.jl) runs the sweep and writes
  the CSV file. It needs ModelingToolkit and takes about 3 minutes.
- [examples/plot_compression.jl](examples/plot_compression.jl) reads that CSV back and does
  everything else: the two checks (`report_sign`, `check_formula`), the analytical formula
  (`analytic_force`) and both figures. It needs nothing but GLMakie, so it loads in a
  second — the plots and the formula can be reworked without re-running the sweep.

`test_compression.jl` includes `plot_compression.jl` at the end, so a full run still
produces the figures and there is only one implementation of each.

`main()` in `test_compression.jl` covers all four points at once: the unstretched lengths 1, 3, 10 and 30 m, the
diameters 2, 4, 6 and 8 mm and the wind speeds 10, 20 and 30 m/s, with
`l_tether_unstretched / l_tether` swept over 0.99 … 1.10, always with 6 segments — 672
operating points, which take about 3 minutes. There is no separate function for point 1 any more: point 2 is a superset
of it, and its CSV carries the same columns, so a dedicated 6 m run added a `mtkcompile`
and a file without adding information.

- `plot_lengths` plots `|mean axial force|` over the relative compression with a
  logarithmic y axis, as a grid of panels — one row per wind speed, one column per
  diameter — with one curve per length. All panels share their axes, so the effect of the
  wind is the shift down the rows and the effect of the diameter the shift across the
  columns.
- `plot_distance` drills into one `(length, diameter, wind)` slice — `slice(results;
  l_unstretched=3.0, d_tether=4.0, v_wind=10.0)` picks one — and plots `|segment force|`
  for every segment and `|anchor force|` for both anchors, also logarithmic. This is what
  showed the six segments carry the same force to within 0.3%.
- The ratio grid (`ratios_around_one`) is geometric in the distance from one, not uniform —
  see below.

  It is the **unstretched** length that is held at 1, 3, 10 and 30 m, and the distance
  `l_tether = l0 / ratio` between the end points that varies, not the other way round. The
  unstretched length is baked into the model and the anchor position is not, so this costs
  one `mtkcompile` per length instead of one per operating point. The swept ratios are
  identical either way and both lengths are in the CSV file, so nothing is lost.

The grid is built from two lists of strain steps, `EXTENSION_STEPS` and
`COMPRESSION_STEPS`. Neither figure uses `MakieControlPlots`: version 0.1.16 has `xscale` but
no `yscale`, so both are built with Makie directly — via `import GLMakie` with every call
qualified, because `using` it as well as `MakieControlPlots` makes their common export
`plot` ambiguous in `Main`, which breaks every example included afterwards.

`main()` writes `data/compression_force_vs_length.csv` with the columns

    v_wind, d_tether, l_tether_unstretched, l_tether, f_top, f_bot, f_mean, f_seg_1 … f_seg_6

It goes to `data` and not to `output`, which is in `.gitignore`: this file is the input of
step two, so it has to survive and be diffable. `data/compression_force.csv` is the
original point-1 run (6 m, 10 m/s) that the table below quotes; nothing regenerates it now.

Open questions and decisions:

- Point 2 originally read "vary l_tether_unstretched between 1.01 and 0.9 of l_tether",
  which literally is 1% compression to 10% extension, i.e. the opposite direction of
  point 1. Resolved in favour of keeping point 1's physical range, and the text of point 2
  corrected accordingly: the ratio `l_tether_unstretched / l_tether` runs from 0.99 (1%
  extension) to 1.10 (10% compression).
- The logarithmic axis shows the magnitude, so that a sign change cannot break it: a curve
  diving towards the clamp `min_force = 1e-4 N` is a force passing through zero.
  `report_sign` prints every operating point whose mean axial force is not tensile, which
  is the direct test of the claim at the top of this document.
- No gravity means the buckled shape is an unstable equilibrium of the spring forces alone,
  so the steady state solver is seeded with a half sine bow in the wind direction whose
  amplitude takes up the slack exactly (`bow_amplitude`).
- Two things in `src/TetherComponent.jl` became parameters, so that a whole strain and wind
  sweep runs on one compiled model: `FixedEnd` holds its node at the position of `pos_fix`
  rather than at a literal, and the wind of `Tether` is `v_wind` rather than
  `se.v_wind_tether`. Only those two, plus the initial states `tether.pos_in` / `vel_in`,
  change between operating points, so the whole run is 5 `mtkcompile` calls for its 182
  operating points. All of them still default to what was passed in, so nothing else
  changes. The tether cross section is parameterised the same way, as `d_tether`,
  `c_spring_unit`, `damping_unit` and `mass_per_m` — the four places the diameter enters —
  so the whole run is 4 `mtkcompile` calls, one per length, for 672 operating points.
- Only the unstretched length still forces a rebuild, because `l_spring(se)` puts `se.l0`
  into the equations as a literal. That is why point 2 holds `l_tether_unstretched` fixed
  and varies the distance, and it is the obvious next parameter if more lengths are wanted.

Note that only two of the four diameter-dependent quantities can move the *equilibrium*:
the stiffness (∝ d²) and the drag area (∝ d). The mass and the damping only shape how the
solver gets there — at rest the accelerations and the spring velocities are zero, so they
drop out of the force balance. They are parameterised anyway, to keep the component
coherent.

### First results (point 1, 6 x 1 m, 10 m/s)

`data/compression_force.csv`, mean axial force per segment, tension positive, with
`r = l_tether_unstretched / l_tether`:

| r | 0.990 | 0.995 | 1.000 | 1.005 | 1.010 | 1.111 |
|---|---|---|---|---|---|---|
| force [N] | 6146 | 3073 | 36.7 | 3.98 | 2.79 | 0.77 |

Two things follow from this.

1. **The force stays tensile over the whole range.** The drag bows a "compressed" tether
   out until the arc is longer than its unstretched length, so the segments are stretched
   even when the end points are closer together than the tether is long. The force never
   changes sign, which is consistent with the claim at the top — the sign it keeps is
   positive, i.e. the tether is never actually in compression.
2. **The spacing had to change.** Below one the stiff tension branch makes the force
   proportional to the distance from one (6146 N at 0.990, 3073 N at 0.995), and it drops
   four decades between 0.990 and 1.005, while above 1.01 it only creeps from 2.8 N to
   0.8 N. A uniform `0.99:0.01:1.10` grid spends 11 of its 12 points on that flat tail and
   resolves none of the interesting part, so the steps shrink towards one instead:
   `EXTENSION_STEPS = [0.01, 0.005, 0.002, 0.001, 0.0005, 0.0002]` and
   `COMPRESSION_STEPS = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1]`, giving 14 points per
   `l_tether`. The force plots are logarithmic for the same reason — on a linear axis the
   single 6146 N point flattens everything else onto the zero line.

### Results of the full sweep (672 points, 3 minutes)

`report_sign` finds **no non-tensile point** anywhere in the sweep. The minimum mean axial
force over all 672 operating points is **0.068 N**, at the case that minimises it in every
respect at once — thinnest tether, weakest wind, shortest length, deepest compression
(2 mm, 10 m/s, `l0 = 1 m`, ratio 1.10). The claim at the top of this document holds across
a 30x range of lengths, 4x of diameters and 3x of wind speeds; the sign it keeps is
positive, i.e. the tether is never actually in compression.

The diameter scaling is remarkably clean. Holding wind, length and ratio fixed and dividing
by the 2 mm value:

| ratio | 4 mm / 2 mm | 6 mm / 2 mm | 8 mm / 2 mm | exponent |
|---|---|---|---|---|
| 0.998 (extension)   | 4.00 | 9.00 | 16.00 | `d²`   |
| 1.000 (zero strain) | 2.52 | 4.33 | 6.35  | `d⁴ᐟ³` |
| ≥ 1.02 (compression)| 2.00 | 3.00 | 4.00  | `d¹`   |

Three regimes, and the exponents are exact to the digits printed (`2^(4/3) = 2.5198`,
`3^(4/3) = 4.3267`, `4^(4/3) = 6.3496`):

- **Extension** is pure spring: the stiffness goes as the cross section, `d²`, and the drag
  is negligible against it.
- **Compression** is pure drag: the force is set by the drag on the bowed tether, whose
  area goes as `d¹`. The stiffness has dropped out entirely — which makes sense, since
  `rel_compression_stiffness` makes the spring 100x softer there.
- **Zero strain** is the crossover, and `d⁴ᐟ³` is the geometric blend of the two.

The exception is the top-right of the range — 30 m/s on a 30 m tether at ratio 0.998 gives
3.62 / 7.94 / 13.97 rather than 4 / 9 / 16 — where the drag is already large enough at 0.2%
extension to pull the exponent away from 2 towards 1.

## Step two: Derive the formula

Done. `analytic_force` in [examples/plot_compression.jl](examples/plot_compression.jl)
reproduces the measured mean axial force of all 672 operating points with a **median error
of 0.003% and a worst case of 0.61%**, with no fitted constant.

### Derivation

A cable under a transverse load that is uniform along its chord hangs in a parabola, and a
chain of straight segments under uniform point loads is exactly that parabola sampled at
its nodes. With `L` the distance between the anchors, `L0` the unstretched length,
`r = L0/L`, `EA = c_spring` the axial stiffness, `n` the number of segments and

    w = ½ ρ c_d d v²        the drag per unit length [N/m]

three relations close the system:

| | |
|---|---|
| sag, from the transverse force balance | `s = w L² / (8F)` |
| arc length of the sampled parabola     | `ΔS = (1 - 1/n²) · 8s²/(3L)` |
| elastic law                            | `ΔS = L0 (1 + F/EA) - L` |

Eliminating `s` and `ΔS` gives a cubic in `F` whose linear term vanishes:

    (r/EA) F³ + (r-1) F² = (1 - 1/n²) · w² L² / 24

It has exactly one positive root — the left side is negative at `F = 0` and grows without
bound — which `analytic_force` returns in closed form (Cardano, with the trigonometric
branch when all three roots are real).

`1 - 1/n²` is the only discretisation term and it is derived, not fitted: an `n`-segment
polyline through a parabola is that much shorter than the smooth curve, so it needs that
much more sag, and hence less force, to take up the same slack. For `n = 6` it is `35/36`.
Dropping it costs a factor of about 1.4% — exactly the residual that was left without it.

### The three regimes fall out of the cubic

Each measured exponent of the diameter is one limit of the same equation:

| regime | limit of the cubic | force | diameter |
|---|---|---|---|
| extension, `r < 1` | drag term negligible | `F = EA (1-r)/r` | `d²` |
| zero strain, `r = 1` | quadratic term vanishes | `F = ((1-1/n²) EA w² L² / 24)^(1/3)` | `d⁴ᐟ³` |
| compression, `r > 1` | elastic term negligible | `F = w L / sqrt(24 (r-1) / (1-1/n²))` | `d¹` |

The middle row is the crossover, and it is where the measured `d⁴ᐟ³` comes from: `EA ∝ d²`
and `w² ∝ d²`, so the cube root of their product goes as `d⁴ᐟ³`. Checked against the data
at `r = 1`, `l0 = 30 m`, 10 m/s, the closed form gives 42.571 / 107.271 / 184.193 /
270.307 N for 2 / 4 / 6 / 8 mm against measured 42.571 / 107.272 / 184.193 / 270.308 N.

The compression branch also explains the rest of the scaling seen in the sweep: `F ∝ v²`
(measured 4.00 and 9.00 going from 10 to 20 and 30 m/s), `F ∝ L`, and `F ∝ 1/sqrt(r-1)`.

### Why the force never changes sign

The compression branch is the answer to the question at the top of this document. Once the
end points are closer together than the tether is long, the drag has to bow the tether out
until its arc is *longer* than its unstretched length — a slack cable cannot carry a
transverse load at all, so it takes up whatever sag the load demands and then some. The
tether is therefore in tension for every `r > 1`, and

    F → w L / sqrt(24 (r-1) / (1-1/n²))

is strictly positive for any non-zero wind. The force only reaches zero in the limit
`w → 0`, i.e. when the tangential wind vanishes — which is the caveat the document opens
with, now quantified.

### Remaining error

The worst 0.61% sits entirely at `r = 1.10`, the deepest compression, where the sag reaches
about 19% of the span and both the parabola and the "drag is uniform along the chord"
assumption start to strain — the model's drag uses the component perpendicular to each
segment, which falls off as the end segments tilt. The residual is one-sided
(`-0.61% … 0.00%`), so it is a systematic shortfall of the approximation, not scatter.
