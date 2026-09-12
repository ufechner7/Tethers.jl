# Julia vs. Python performance

Both ecosystems can give a stiff tether model an analytic, symbolically
generated, sparse Jacobian, and neither is fast without one. This page compares
them on the `Tether_08` model (ten segments unless stated otherwise, tether
drag, reel-out, one fixed and one free end point), over 10 s, sampled every 20
ms, at a relative and absolute tolerance of $10^{-6}$, from a straight-line
initial condition. CPU: Ryzen 9 7950X, Julia 1.13, ModelingToolkit 11, CasADi
3.8.

|            | Julia                          | Python                                           |
| ---------- | ------------------------------ | ------------------------------------------------ |
| Model      | ModelingToolkit `System`       | CasADi `SX` expression graph                     |
| Jacobian   | `ODEProblem(sys, …; jac=true)` | `ca.jacobian(ydot, y)`                           |
| Sparsity   | `sparse=true`                  | found by CasADi, always                          |
| Integrator | `FBDF`                         | SUNDIALS CVODES, `linear_multistep_method='bdf'` |

`FBDF` is a fixed-leading-coefficient BDF, which is the variant SUNDIALS
implements, so CVODES-BDF is the like-for-like match. IDAS solves the same
formulas for an implicit DAE; it is 3-11% slower than CVODES here, so the
examples use CVODES.

## The Jacobian dominates

Solve time in ms, median and interquartile range over the samples
`BenchmarkTools` fits into twelve seconds per configuration:

| segments | states | Julia AD, dense |   Julia `jac` | Julia `jac+sparse` | CasADi CVODES sparse |
| -------: | -----: | --------------: | ------------: | -----------------: | -------------------: |
|        5 |     36 |       13.0 ±0.4 |     12.7 ±0.7 |           9.1 ±1.0 |            12.3 ±0.5 |
|       10 |     66 |       31.9 ±1.1 |     27.9 ±0.2 |          17.1 ±1.9 |            20.7 ±0.2 |
|       20 |    126 |      485.0 ±6.1 |    438.4 ±6.1 |        177.9 ±31.6 |           170.0 ±4.9 |
|       40 |    246 |    4332.9 ±22.5 | 3765.9 ±154.5 |        603.1 ±10.5 |          617.1 ±13.8 |

CasADi with a dense linear solver instead costs 13.9, 28.7, 1037.2 and 4551.8
ms, so the sparse factorization is worth 7.4x at forty segments there. The Julia
column comes from `ODEProblem(sys, …; jac, sparse)`, the CasADi one from
`examples/python/bench_casadi.py`.

The Jacobian is block-tridiagonal: 591 of 4356 entries at ten segments, 2391 of
60516 at forty. The analytic Jacobian alone buys little, because the dense
factorization then dominates; it is the two together that pay, by 1.4x at five
segments and 7.2x at forty; the same holds on the CasADi side.

With `jac=true` the `autodiff` keyword no longer affects anything: the solver
uses the supplied Jacobian and never differentiates. `FBDF()` and
`FBDF(autodiff=AutoForwardDiff())` measure the same within noise.

Head to head, with both sides analytic and sparse, Julia's `FBDF` and SUNDIALS
CVODES-BDF are within 20% of each other, Julia ahead on the small models and
CVODES on the large ones. CVODES beats IDAS at every size here (10.8 against
13.2 ms at five segments, 601 against 703 at forty), so the implicit-DAE
formulation costs a little.

Generating the Jacobian is not free. ModelingToolkit needs 0.5 s at five
segments and 4.1 s at forty, on top of `mtkcompile`; `ca.jacobian` needs 2 ms
and 16 ms. For a script that solves once, that build cost outweighs the saving;
it pays back when a compiled model is re-solved.

## Where the old numbers came from

Earlier versions of this page reported Julia as 13 to 30 times faster than
Python. That compared compiled Julia against a hand-derived Jacobian evaluated
in interpreted NumPy, which costs 707 µs per call against 33 µs for the same
Jacobian as a CasADi function — a factor of 21 in the inner loop of every Newton
iteration. It measured the binding, not the language. With both sides compiled
and sparse, they are within 20% of each other.

The code size changed with it. The Python examples used to carry the Jacobian of
every model as hand-written calculus - 86 to 150 lines per file in
`Tether_06.py`, `Tether_06c.py`, `Tether_07.py` and `Tether_08.py`.
`ca.jacobian` replaced all of it, and the ten examples lost about 700 lines
between them. The Julia examples never wrote a Jacobian at all.

## Caveats

- `FBDF` and CVODES are different BDF implementations, with their own step-size
  and order heuristics; they take different paths through the same problem.
- The CasADi timings cross into Python once per solve, not once per Newton
  iteration, because the whole model is compiled into the integrator. A model
  driven from a Python callback per step pays an interpreter cost that the table
  above avoids.
- The tether's taut/slack switch is a hard step. Its exact derivative is zero on
  either side, which is what ModelingToolkit, CasADi's `if_else` and ForwardDiff
  all report. Finite differences instead return a large secant slope across the
  step, and some configurations only converge because of it — see the
  steady-state solves of `Tether_08`–`Tether_11`, and `Tether_11`'s time
  simulation.
