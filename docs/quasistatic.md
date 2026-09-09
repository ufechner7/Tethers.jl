# Quasi-static tether model: state of the port

Notes on porting the quasi-static tether model from the `andrea_quasistatic`
branch onto `main`, and on the open questions that port uncovered.

## Why the branch was replanted instead of rebased

`andrea_quasistatic` forked 115 commits before `main` was restructured. In the
meantime `main` moved the examples from `src/` to `examples/`, upgraded to
ModelingToolkit 11, replaced ControlPlots with MakieControlPlots and moved to
Julia 1.11/1.12. Of the 85 distinct commit subjects on the branch, 75 do not
exist on `main`, and 86 of the commits are old-layout work touching Manifest
files and example paths that `main` has since deleted or moved.

Replaying that history would have meant re-resolving the same layout conflicts
on almost every commit. The feature itself is small, so it was replanted onto
`main` as a fresh set of commits instead. The resulting diff is additive:
roughly 1500 inserted lines and one deleted, against the original pull
request's 4161 deletions of Manifest churn.

### What moved where

| Branch | Now | Note |
| --- | --- | --- |
| `src/Tether_quasistatic.jl` | unchanged | numerics untouched |
| `src/Tether_qsm_dual.jl` | unchanged | numerics untouched |
| `examples_quasistatic/` | `examples/quasistatic/` | shares `examples/Project.toml` |
| `src/Tether_10.jl` | `examples/Tether_11.jl` | renamed, see below |
| `test/test_qsm.jl`, `test/data/` | unchanged | wired into `runtests.jl` |

`src/Tether_10.jl` had to be renamed because `main` already has an unrelated
`examples/Tether_10.jl` (the re-usable acausal component). The branch's file is
the imposed-kite-motion example, so it became `examples/Tether_11.jl` and was
added to `examples/menu.jl`.

The six example scripts used the Matplotlib-style `plt.` API that ControlPlots
exposed. MakieControlPlots has no such passthrough, so they were rewritten
against GLMakie directly, following the `import GLMakie` idiom already used by
`examples/Tether_09.jl`. All six run.

## Open question: the angle convention

**This is the one that needs a decision.**

`res!` in `src/Tether_quasistatic.jl` and the MATLAB reference data in
`test/data/` disagree on how the state vector's two angles define the tether
direction at the ground station.

`res!` uses an elevation/azimuth convention, where `θ` is measured up from the
x–y plane:

```julia
dir ∝ [cos(θ)cos(φ), cos(θ)sin(φ), sin(θ)]
```

The reference data was produced with a z-up convention, where `θ` is measured
from the vertical:

```julia
dir ∝ [sin(θ)cos(φ), sin(φ), cos(θ)cos(φ)]
```

For `test/data/input_basic_test.mat` (`θ = 18.435°`, `φ = -17.548°`,
`Tn = 160941 N`, kite at `[100, 100, 300]`, 15 segments) the two give first
segments pointing in quite different directions:

| | first segment direction | `p0` |
| --- | --- | --- |
| `res!` as written | `[0.905, -0.286, 0.316]` | `[391.18, -123.70, 136.76]` |
| reference data | `[0.302, -0.302, 0.905]` | `[129.36, -129.36, 391.88]` |

The computed tether is a differently-oriented one, not a slightly inaccurate
one: `‖p0 - p0_ref‖ = 365.6 m` on a 431 m tether.

Re-running `res!` with the angles converted into the reference convention drops
that to **1.61 m**, which is what identifies the convention as the cause. A
third candidate, the standard spherical form `[sinθcosφ, sinθsinφ, cosθ]`, is
ruled out at 90.0 m. Note that `[sin(θ)cos(φ), sin(φ), cos(θ)cos(φ)]` and
`[sin(θ), tan(φ), cos(θ)]` are the same vector once normalised, and give
identical results.

### Tolerance will not paper over this

| | max relative error | rtol needed to pass |
| --- | --- | --- |
| as written | 202 % | `2.02` |
| angles converted | 0.8 % | `0.008` |

At `rtol = 2.02` the assertion would no longer constrain anything. Only after
the convention is settled does a tolerance become meaningful — and then a
modest one does the job, though the residual 1.6 m is a second, smaller
discrepancy that still wants explaining.

### Two defensible resolutions

1. `res!` has the wrong convention and should be changed to match the
   reference. This alters what `state_vec` means for every caller, so
   `init_quasistatic`, `simulate_tether` and all six examples are affected.
2. The `.mat` files are simply expressed in the MATLAB frame, and
   `get_initial_conditions` should convert on load, leaving `res!` alone.

Which is correct depends on which frame the model is meant to expose. It was
left alone pending that decision; the four value comparisons in
`test/test_qsm.jl` are marked `@test_broken` with the diagnosis in a comment,
so the suite is green (4 pass, 4 broken) and the reference data keeps its
purpose.

### Reproducing

```julia
include("src/Tether_quasistatic.jl")
sv, kp, kv, wv, tl, se = get_initial_conditions("test/data/input_basic_test.mat")
ref = matread("test/data/basic_test_results.mat")
Ns = size(wv, 2)
buffers() = [zeros(3, Ns) for _ in 1:5]

_, _, _, p0 = res!(zeros(3), sv, (kp, kv, wv, tl, se, buffers(), Ns, true))
norm(p0 .- vec(ref["p0"]))          # 365.6

# the same call with the angles converted to the reference convention
d = [sin(sv[1])cos(sv[2]), sin(sv[2]), cos(sv[1])cos(sv[2])]; d ./= norm(d)
sv2 = MVector(asin(d[3]), atan(d[2], d[1]), sv[3])
_, _, _, p0c = res!(zeros(3), sv2, (kp, kv, wv, tl, se, buffers(), Ns, true))
norm(p0c .- vec(ref["p0"]))         # 1.61
```

## Bugs found and fixed

**The test had never run.** `test/test_qsm.jl` passed `res!` a seven-element
`param` tuple, but `res!` destructures eight values from it, so the call died
with a `BoundsError` before reaching any assertion. This is why the convention
mismatch above went unnoticed for so long — it is a pre-existing discrepancy,
not a regression from the port. The test now builds the eight-field named
tuple, and locates its `.mat` files relative to `@__DIR__` rather than the
current working directory.

It also asserted `p0 isa Vector`, which is false: `res!` returns an
`MVector{3,Float64}`. That assertion was corrected rather than deleted.

**`Tether_11.jl` was entirely dead code.** Its `main()` and the call that drove
it sat inside a `"""..."""` string literal, so including the file defined
`model` and `simulate` and then did nothing at all. `main()` was lifted out and
the plotting block ported to GLMakie.

**`Tether_11.jl` did not build under MTK 11.** Writing the differential
equations as

```julia
eqs1 = vcat(D.(pos) .~ vel, D.(vel) .~ acc)
eqs2 = vcat(eqs1...)
```

leaves a `Vector{Any}`, which the ModelingToolkit 11 `System` constructor
rejects. It now builds them per column, the way `examples/Tether_08.jl` does.

**A stale compat bound.** `PreallocationTools = "0.4.25"` was carried over from
the branch and predates that package's 1.0 release, so resolving the workspace
failed with an empty intersection against 1.1.2. Widened to `"0.4.25, 1"`. The
bound only bites from the workspace root, where every member project's compat
is reconciled at once — resolving `--project=test` alone did not surface it.

## Open question: the `maxiters` warning

`examples/Tether_11.jl` emits this during its steady-state solve, while
`examples/Tether_08.jl` does not:

```
Interrupted. Larger maxiters is needed. ...
```

The example completes and plots correctly regardless; the warning comes from
the `DynamicSS` call that produces the initial tether shape.

The obvious hypothesis is that `model` passes the prescribed `acc_p2` straight
into the steady-state solve, so the second end point accelerates forever and no
steady state exists — the same reason `se.v_ro` is zeroed on the line above.
**This was tested and is wrong**: zeroing `acc_p2` for the steady-state solve
does not remove the warning. The change was reverted, and the real cause is
still unknown.

## Verified

- The workspace root and the `test` environment both resolve against MTK 11.
- `test/test_qsm.jl`: 4 pass, 4 broken.
- All six `examples/quasistatic/` scripts run.
- `examples/Tether_11.jl` runs end to end (~67 s).
