# -*- coding: utf-8 -*-
"""
Benchmark of the Tether_08 model with a CasADi-generated Jacobian, for comparison with
the ModelingToolkit version of examples/Tether_08.jl.

Tether_08.py derives the Jacobian of the particle accelerations by hand, in about a
hundred lines of calculus (calc_segment_force_jac, calc_acc_jac). CasADi derives the same
Jacobian from the model expression itself:

    J = ca.jacobian(ydot, y)

This script checks that the two agree, reports the sparsity CasADi finds, and times the
solve with the SUNDIALS integrators that CasADi ships, so that both sides of the
Julia/Python comparison use an analytic, symbolically generated, sparse Jacobian.

CVODES with linear_multistep_method='bdf' is the like-for-like match for the FBDF solver
of the Julia examples: both are fixed-leading-coefficient BDF with a Newton iteration.
IDAS solves the same formulas for an implicit DAE and is what Tether_08.py uses.

Not part of the tutorial and not needed to run the examples. Requires casadi, which is
not in CondaPkg.toml:

    pip install casadi
"""
import math
import statistics
import timeit

import numpy as np
import casadi as ca

from dataclasses import dataclass, field


@dataclass
class Settings3:
    """Same physics and same units as Settings3 of Tether_08.py."""
    g_earth: np.ndarray = field(default_factory=lambda: np.array([0.0, 0.0, -9.81]))
    v_wind_tether: np.ndarray = field(default_factory=lambda: np.array([2.0, 0.0, 0.0]))
    rho: float = 1.225
    cd_tether: float = 0.958
    l0: float = 70.0
    v_ro: float = 0.3
    d_tether: float = 4.0
    rho_tether: float = 724.0
    c_spring: float = 614600.0
    rel_compression_stiffness: float = 0.01
    damping: float = 473.0
    segments: int = 10
    duration: float = 10.0


def build(se):
    """ The whole model as one CasADi expression graph, mirroring calc_accelerations of
        Tether_08.py: p1 is held fixed, p2 is free. Returns the symbolic time, the state
        vector `y` = (pos, vel), its derivative `ydot`, and the straight-line initial
        state between the two end points. """
    n = se.segments
    split = 3 * (n + 1)                         # start of the velocities in y
    t_sym = ca.SX.sym('t')
    y = ca.SX.sym('y', 2 * split)
    pos = [y[3*i:3*i + 3] for i in range(n + 1)]
    vel = [y[split + 3*i:split + 3*i + 3] for i in range(n + 1)]

    l_spring = (se.l0 + se.v_ro * t_sym) / n
    c_spring = se.c_spring / l_spring
    damping = se.damping / l_spring
    m_particle = se.rho_tether * math.pi * (se.d_tether / 2000.0) ** 2 * l_spring

    force = [ca.SX.zeros(3) for _ in range(n + 1)]
    for j in range(n):
        segment = pos[j + 1] - pos[j]
        length = ca.norm_2(segment)
        e = segment / length
        spring_vel = ca.dot(e, vel[j + 1] - vel[j])
        rcs = se.rel_compression_stiffness
        # the same hard taut/slack switch as Tether_08.jl; if_else is differentiated on
        # the active branch, which is the convention Tether_08.py's hand-derived
        # Jacobian uses as well
        c_spr = c_spring / (1.0 + rcs) * (rcs + ca.if_else(length > l_spring, 1.0, 0.0))
        fs = (c_spr * (length - l_spring) + damping * spring_vel) * e

        v_app = se.v_wind_tether - (vel[j] + vel[j + 1]) / 2.0
        v_app_perp = v_app - ca.dot(v_app, e) * e
        drag = (0.25 * se.rho * se.cd_tether * ca.norm_2(v_app_perp)
                * length * se.d_tether / 1000.0) * v_app_perp
        force[j] = force[j] + fs + drag
        force[j + 1] = force[j + 1] - fs + drag

    acc = []
    for i in range(n + 1):
        if i == 0:
            acc.append(ca.SX.zeros(3))          # p1 is fixed
        else:
            mass = 0.5 * m_particle if i == n else m_particle
            acc.append(se.g_earth + force[i] / mass)

    p1 = np.array([0.0, 0.0, 0.0])
    p2 = np.array([-40.0, 0.0, -47.0])
    pos0 = np.array([p1 + i * (p2 - p1) / n for i in range(n + 1)])
    y0 = np.concatenate([pos0.flatten(), np.zeros(split)])
    return t_sym, y, ca.vertcat(*vel, *acc), y0


# CVODES-BDF is the match for the Julia examples' FBDF; IDAS is what Tether_08.py uses.
# The dense variants are the analogue of ODEProblem(..., jac=true) without sparse=true.
def sample(call, seconds=8.0):
    """ Median, interquartile range and sample count of `call`, in ms. One call per sample,
        repeated until `seconds` of wall clock have been spent, so that a slow
        configuration is not sampled to death and a fast one still gets many samples. """
    call()                                                  # warm up
    once = min(timeit.repeat(call, number=1, repeat=3))
    reps = max(7, min(500, int(seconds / max(once, 1e-6))))
    times = sorted(t * 1000 for t in timeit.repeat(call, number=1, repeat=reps))
    q1, q3 = statistics.quantiles(times, n=4)[0], statistics.quantiles(times, n=4)[2]
    return statistics.median(times), q3 - q1, len(times)


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
DT = 0.02


def check_against_handwritten(se, t_sym, y, ydot, y0):
    """ Compare the CasADi Jacobian with the hand-derived calc_acc_jac of Tether_08.py at
        a generic state. Skipped if Tether_08.py cannot be imported (it needs Assimulo). """
    try:
        from Tether_08 import Settings3 as S8, set_tether_diameter, calc_accelerations, calc_acc_jac
    except Exception as exc:
        print(f"  (skipped, cannot import Tether_08.py: {type(exc).__name__}: {exc})")
        return
    n = se.segments
    split = 3 * (n + 1)
    f_rhs = ca.Function('f', [t_sym, y], [ydot])
    f_jac = ca.Function('J', [t_sym, y], [ca.jacobian(ydot, y)])

    rng = np.random.default_rng(0)
    y_test = y0 + np.concatenate([rng.normal(0, 0.5, split), rng.normal(0, 1.0, split)])
    t_test = 3.7
    pos = y_test[:split].reshape(n + 1, 3)
    vel = y_test[split:].reshape(n + 1, 3)

    se8 = S8()
    set_tether_diameter(se8, se.d_tether)
    acc_ref = calc_accelerations(t_test, pos, vel, se8, True, False, se8.v_ro)
    dacc_dpos, dacc_dvel = calc_acc_jac(t_test, pos, vel, se8, True, False, se8.v_ro)
    acc_cas = np.array(f_rhs(t_test, y_test)).flatten()[split:].reshape(n + 1, 3)
    jac_cas = np.array(f_jac(t_test, y_test))

    print(f"  max |acc         - hand-derived| = {np.max(np.abs(acc_cas - acc_ref)):.3e}")
    print(f"  max |d(acc)/dpos - hand-derived| = {np.max(np.abs(jac_cas[split:, :split] - dacc_dpos)):.3e}")
    print(f"  max |d(acc)/dvel - hand-derived| = {np.max(np.abs(jac_cas[split:, split:] - dacc_dvel)):.3e}")


def main(segment_counts=(5, 10, 20, 40)):
    se = Settings3()
    print(f"Jacobian cross-check against Tether_08.py, {se.segments} segments:")
    t_sym, y, ydot, y0 = build(se)
    check_against_handwritten(se, t_sym, y, ydot, y0)

    print(f"\n{'n':>3} {'states':>7} {'jac nnz':>8} {'%dense':>7} {'jac gen':>9}  "
          + "".join(f"{name:>28}" for name, _, _ in CONFIGS))
    print(f"{'':>3} {'':>7} {'':>8} {'':>7} {'':>9}  "
          + "".join(f"{'median ±IQR [ms], n':>28}" for _ in CONFIGS))
    for n in segment_counts:
        se = Settings3(segments=n)
        t_sym, y, ydot, y0 = build(se)
        t_jac = min(timeit.repeat(lambda: ca.jacobian(ydot, y), number=1, repeat=5))
        jac = ca.jacobian(ydot, y)

        grid = np.arange(0.0, se.duration + DT / 2, DT)
        dae = {'x': y, 't': t_sym, 'ode': ydot}
        row = []
        for _, plugin, opts in CONFIGS:
            options = dict(opts, abstol=1e-6, reltol=1e-6)
            integrator = ca.integrator('I', plugin, dae, 0.0, grid[1:], options)
            med, iqr, reps = sample(lambda: integrator(x0=y0))
            row.append(f"{med:>15.1f} ±{iqr:<6.1f} n={reps:<4d}")
        size = 2 * 3 * (n + 1)
        print(f"{n:>3} {size:>7} {jac.sparsity().nnz():>8} "
              f"{100 * jac.sparsity().nnz() / size ** 2:>6.1f}% {t_jac * 1000:>8.1f}ms  "
              + "".join(row))


if __name__ == '__main__':
    main()
