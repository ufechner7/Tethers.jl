# Julia vs. Python performance — Tether_06

Comparison of the two `Tether_06` implementations, both solving the same DAE (5-segment
damped tether, taut/slack spring with a smoothed 1% compression floor) with an exact
analytic Jacobian.

| | Julia (`FBDF`, ModelingToolkit) | Python (`IDA`, Assimulo) |
|---|---|---|
| Solve time | ~0.002–0.003s | ~0.041s |
| Output points | 501 (`saveat`, exact) | 501 (resampled onto exact grid) |
| Steps taken | not reported by SciML by default | 436 |
| Jacobian evals | analytic, symbolically generated (compiled once, ~free per call) | 83 (analytic, hand-derived) |
| Setup/compile cost | `mtkcompile` symbolic pipeline (one-time, amortized across repeat solves — not counted here) | negligible (Python interpreted, no compile step) |

**Julia is about 15–20x faster** than Python for this problem, even with both now using
exact analytic Jacobians.

## Why the gap remains

Despite both being "Newton + analytic Jacobian" DAE solvers:

- **Compiled vs. interpreted inner loop**: ModelingToolkit generates and compiles native
  Julia functions for the residual and Jacobian ahead of time; Python's `res`/`jac`
  methods re-run interpreted NumPy code (with per-call array allocation, Python-level
  loops, function-call overhead) on every single Newton iteration.
- **Allocation overhead**: Julia's solve reused ~5,000 small allocations total across the
  whole 10s run; Python's `res`/`jac` allocate fresh NumPy arrays (segment vectors, outer
  products, etc.) on every call — 436 steps × several evaluations per step × several
  small array allocations each adds up in Python's interpreter/GC overhead in a way it
  doesn't in compiled Julia.
- **Fixed cost per Python↔C boundary crossing**: Assimulo's IDA is itself compiled
  Sundials C code, but every residual/Jacobian evaluation has to cross back into the
  Python interpreter, whereas Julia's whole pipeline (solver + residual + Jacobian) is
  native compiled code with no language-boundary crossings.

## Context

Both are fast in absolute terms for this problem (2ms vs 41ms is imperceptible either
way). The practically important result is that the Python side went from **79.2s to
0.04s** (~2000x) over the course of this optimization work, closing nearly all of the
gap that existed purely from using a finite-difference Jacobian instead of an analytic
one.
