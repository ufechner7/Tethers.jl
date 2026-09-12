# -*- coding: utf-8 -*-
"""
Benchmark of the Tether_08 model, for the comparison in docs/julia_vs_python.md.

The model itself comes from Tether_08.py - this script only varies the number of segments
and the integrator, so that there is one definition of the physics and not two. It reports
how the Jacobian's sparsity grows with the segment count, what CasADi spends deriving it,
and how the SUNDIALS integrators compare with a sparse and with a dense linear solver.

CVODES with linear_multistep_method='bdf' is the like-for-like match for the FBDF solver of
the Julia examples: both are fixed-leading-coefficient BDF with a Newton iteration. IDAS
solves the same formulas for an implicit DAE.

Not part of the tutorial; run it directly:

    python examples/python/bench_casadi.py
"""
import statistics
import timeit

import numpy as np
import casadi as ca

from Tether_08 import Settings3, build_model, calc_initial_state, set_tether_diameter

# the endpoints Tether_08.py and Tether_08.jl both use
P1 = np.array([0.0, 0.0, 0.0])
P2 = np.array([-40.0, 0.0, -47.0])
DT = 0.02
DURATION = 10.0

CONFIGS = [
    ("cvodes BDF, sparse LU", 'cvodes', {'linear_multistep_method': 'bdf',
                                         'nonlinear_solver_iteration': 'newton',
                                         'linear_solver': 'csparse'}),
    ("cvodes BDF, dense LU",  'cvodes', {'linear_multistep_method': 'bdf',
                                         'nonlinear_solver_iteration': 'newton',
                                         'linear_solver': 'lapacklu'}),
    ("idas, sparse LU",       'idas',   {'linear_solver': 'csparse'}),
    ("idas, dense LU",        'idas',   {'linear_solver': 'lapacklu'}),
]


def sample(call, seconds=8.0):
    """ Median, interquartile range and sample count of `call`, in ms. One call per sample,
        repeated until `seconds` of wall clock have been spent, so that a slow
        configuration is not sampled to death and a fast one still gets many samples. """
    call()                                                  # warm up
    once = min(timeit.repeat(call, number=1, repeat=3))
    reps = max(7, min(500, int(seconds / max(once, 1e-6))))
    times = sorted(t * 1000 for t in timeit.repeat(call, number=1, repeat=reps))
    quartiles = statistics.quantiles(times, n=4)
    return statistics.median(times), quartiles[2] - quartiles[0], len(times)


def main(segment_counts=(5, 10, 20, 40)):
    print(f"{'n':>3} {'states':>7} {'jac nnz':>8} {'%dense':>7} {'jac gen':>9}  "
          + "".join(f"{name:>28}" for name, _, _ in CONFIGS))
    print(f"{'':>3} {'':>7} {'':>8} {'':>7} {'':>9}  "
          + "".join(f"{'median ±IQR [ms], n':>28}" for _ in CONFIGS))
    for n in segment_counts:
        se = Settings3(segments=n, duration=DURATION)
        set_tether_diameter(se, se.d_tether)
        pos0, vel0, _, _ = calc_initial_state(se, P1, P2)
        t, y, ydot = build_model(se, fix_p1=True, fix_p2=False, v_ro=se.v_ro)
        y0 = np.concatenate([pos0.flatten(), vel0.flatten()])

        t_jac = min(timeit.repeat(lambda: ca.jacobian(ydot, y), number=1, repeat=5))
        jac = ca.jacobian(ydot, y)

        grid = np.arange(0.0, se.duration + DT / 2, DT)
        dae = {'x': y, 't': t, 'ode': ydot}
        row = []
        for _, plugin, opts in CONFIGS:
            integrator = ca.integrator('bench', plugin, dae, 0.0, grid[1:],
                                       dict(opts, abstol=1e-6, reltol=1e-6))
            med, iqr, reps = sample(lambda: integrator(x0=y0))
            row.append(f"{med:>15.1f} ±{iqr:<6.1f} n={reps:<4d}")
        size = 2 * 3 * (n + 1)
        print(f"{n:>3} {size:>7} {jac.sparsity().nnz():>8} "
              f"{100 * jac.sparsity().nnz() / size ** 2:>6.1f}% {t_jac * 1000:>8.1f}ms  "
              + "".join(row))


if __name__ == '__main__':
    main()
