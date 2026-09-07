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

### Status

Both points are implemented in [examples/test_compression.jl](examples/test_compression.jl)
and the script is reachable from `examples/menu.jl`.

`main()` covers all three points at once: the unstretched lengths 1, 3, 10 and 30 m and the
wind speeds 10, 20 and 30 m/s, with `l_tether_unstretched / l_tether` swept over
0.99 … 1.10, always with 6 segments — 168 operating points. There is no separate function
for point 1 any more: point 2 is a superset of it, and its CSV carries the same columns, so
a dedicated 6 m run added a `mtkcompile` and a file without adding information.

- `plot_lengths` plots `|mean axial force|` over the relative compression with a
  logarithmic y axis, one panel per wind speed and one curve per length, all panels sharing
  their axes so the wind can be read off by comparing them.
- `plot_distance` drills into one `(length, wind)` slice and plots `|segment force|` for
  every segment and `|anchor force|` for both anchors, also logarithmic. This is what
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

    v_wind, l_tether_unstretched, l_tether, f_top, f_bot, f_mean, f_seg_1 … f_seg_6

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
  operating points. Both parameters still default to what was passed in, so nothing else
  changes. With point 1 folded into point 2 this is 4 calls for 168 points.
- Only the unstretched length still forces a rebuild, because `l_spring(se)` puts `se.l0`
  into the equations as a literal. That is why point 2 holds `l_tether_unstretched` fixed
  and varies the distance, and it is the obvious next parameter if more lengths are wanted.

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

## Step two: Derive the formula

Open. Use the CSV files above. The formula has to reproduce two regimes and the crossover
between them: the linear tension branch below `r = 1`, and the drag-dominated branch above
it, where the force is set by the balance between the aerodynamic drag on the bowed tether
and the axial force needed to hold that bow.
