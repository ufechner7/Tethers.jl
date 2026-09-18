# -*- coding: utf-8 -*-
"""
Tutorial example simulating a 3D mass-spring system with a nonlinear spring (1% stiffness
for l < l_0), n tether segments, tether drag and reel-in and reel-out.

The model is written once as a CasADi expression graph. CasADi differentiates it to get
the exact Jacobian and finds its sparsity; the same expression serves the steady-state
solve and the time integration, so the physics has a single definition. SUNDIALS' CVODES
integrates it with a BDF formula and a sparse Newton solve - the same
fixed-leading-coefficient BDF that the FBDF solver of Tether_08.jl uses.

Like Tether_08.jl, the initial tether shape for a given pair of endpoints is found first,
here with scipy.optimize.least_squares, by solving for the positions at which the
acceleration of every free particle is zero (with velocities zero and v_ro = 0).

State vector y = pos[0..n], vel[0..n] (each particle contributes a 3D vector).
"""
import math
import os
import time as time_module
from dataclasses import dataclass, field

import numpy as np
import casadi as ca
import matplotlib.pyplot as plt
from scipy.optimize import least_squares


@dataclass
class Settings3:
    g_earth: np.ndarray = field(default_factory=lambda: np.array([0.0, 0.0, -9.81]))  # [m/s^2]
    v_wind_tether: np.ndarray = field(default_factory=lambda: np.array([2.0, 0.0, 0.0]))
    rho: float = 1.225
    cd_tether: float = 0.958
    l0: float = 70.0                             # initial tether length             [m]
    v_ro: float = 0.3                            # reel-out speed                  [m/s]
    d_tether: float = 4.0                        # tether diameter                  [mm]
    rho_tether: float = 724.0                    # density of Dyneema            [kg/m^3]
    c_spring: float = 614600.0                   # unit spring constant              [N]
    rel_compression_stiffness: float = 0.01      # relative compression stiffness    [-]
    damping: float = 473.0                       # unit damping constant            [Ns]
    segments: int = 10                           # number of tether segments         [-]
    alpha0: float = math.pi / 10                 # initial tether angle            [rad]
    duration: float = 30.0                       # duration of the simulation        [s]
    save: bool = False                           # save png files in folder video


def set_tether_diameter(se, d, c_spring_4mm=614600.0, damping_4mm=473.0):
    """ Set the tether diameter `d` [mm] and scale the unit spring constant and unit
        damping constant with the cross section, relative to a 4 mm reference tether. """
    se.d_tether = d
    se.c_spring = c_spring_4mm * (d / 4.0) ** 2
    se.damping = damping_4mm * (d / 4.0) ** 2


def calc_initial_state(se, p1, p2):
    """ Linearly interpolated initial position (and zero velocity) for each tether
        particle between the endpoints p1 and p2. If one endpoint is None, it is derived
        from the other one using se.alpha0 and se.l0; at least one must be given.
        Returns (pos0, vel0, p1, p2), each pos/vel a (segments+1) x 3 array. """
    if p1 is None and p2 is None:
        raise ValueError("at least one of p1 and p2 must be defined")
    if p2 is None:
        z = math.cos(se.alpha0) * se.l0
        y = math.sin(se.alpha0) * se.l0
        p2 = np.array([p1[0], p1[1] - y, p1[2] - z])
        print("p2:", p2)
    elif p1 is None:
        z = math.cos(se.alpha0) * se.l0
        y = math.sin(se.alpha0) * se.l0
        p1 = np.array([p2[0], p2[1] + y, p2[2] + z])
        print("p1:", p1)
    p1 = np.asarray(p1, dtype=float)
    p2 = np.asarray(p2, dtype=float)

    n = se.segments
    pos0 = np.zeros((n + 1, 3))
    vel0 = np.zeros((n + 1, 3))
    delta = (p2 - p1) / n
    for i in range(n + 1):
        pos0[i, :] = p1 + i * delta
    return pos0, vel0, p1, p2


def add_initial_sag(se, pos0):
    """ Bow the interior particles of the straight-line initial guess `pos0` downwards,
        such that the length of the resulting polyline matches the unstretched tether
        length se.l0. The endpoints are not moved (the bow is zero at both ends).

        This matters because the tether is slack whenever the distance between the two
        endpoints is smaller than se.l0: the straight line is then a uniformly compressed
        chain, and the steady-state solver converges from it to the nearest root, which is
        the *unstable* arch-shaped equilibrium (every segment compressed, bowing upwards)
        instead of the hanging one. Tether_08.jl does not need this, because its
        DynamicSS solver integrates the damped equations of motion forward in time and
        therefore always relaxes into the stable, hanging equilibrium. """
    n = se.segments
    chord = np.linalg.norm(pos0[n] - pos0[0])
    if chord >= se.l0:
        return pos0                     # the tether is taut, the straight line is fine
    # parabolic bow: zero at both endpoints, maximum in the middle
    xi = np.linspace(0.0, 1.0, n + 1)
    shape = 4.0 * xi * (1.0 - xi)

    def polyline_length(sag):
        pos = pos0.copy()
        pos[:, 2] -= sag * shape
        return np.sum(np.linalg.norm(pos[1:] - pos[:-1], axis=1))

    # the polyline length grows monotonically with the sag, so bisection finds the
    # amplitude at which it equals se.l0
    lo, hi = 0.0, se.l0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if polyline_length(mid) < se.l0:
            lo = mid
        else:
            hi = mid
    pos = pos0.copy()
    pos[:, 2] -= 0.5 * (lo + hi) * shape
    return pos


def accelerations(t, pos, vel, se, fix_p1, fix_p2, v_ro):
    """ Acceleration of every tether particle, as a list of CasADi 3-vectors, given the
        positions `pos` and velocities `vel` (lists of 3-vectors), the reel-out speed
        `v_ro` and whether the two endpoints are held fixed. Mirrors the equations built
        by Tether_08.jl's `model` function.

        The spring is nonlinear: a taut segment (length > l_spring) uses the full
        stiffness, a slack one only se.rel_compression_stiffness of it. `if_else` is
        differentiated on the active branch, so the Jacobian is the one of the branch the
        segment is currently on. Only the component of the apparent wind perpendicular to
        a segment produces drag, and half of a segment's drag acts on each of its ends. """
    n = se.segments
    l_spring = (se.l0 + v_ro * t) / n
    c_spring = se.c_spring / l_spring
    damping = se.damping / l_spring
    m_particle = se.rho_tether * math.pi * (se.d_tether / 2000.0) ** 2 * l_spring

    force = [ca.SX.zeros(3) for _ in range(n + 1)]
    for j in range(n):
        segment = pos[j + 1] - pos[j]
        length = ca.norm_2(segment)
        e = segment / length                              # unit vector, j towards j+1
        spring_vel = ca.dot(e, vel[j + 1] - vel[j])       # rate of change of the length
        rcs = se.rel_compression_stiffness
        c_spr = c_spring / (1.0 + rcs) * (rcs + ca.if_else(length > l_spring, 1.0, 0.0))
        fs = (c_spr * (length - l_spring) + damping * spring_vel) * e

        v_app = se.v_wind_tether - (vel[j] + vel[j + 1]) / 2.0
        v_app_perp = v_app - ca.dot(v_app, e) * e
        drag = (0.25 * se.rho * se.cd_tether * ca.norm_2(v_app_perp)
                * length * se.d_tether / 1000.0) * v_app_perp
        force[j] = force[j] + fs + drag       # segment j pulls particle j towards j+1 ...
        force[j + 1] = force[j + 1] - fs + drag   # ... and particle j+1 towards particle j

    acc = []
    for i in range(n + 1):
        if (i == 0 and fix_p1) or (i == n and fix_p2):
            acc.append(ca.SX.zeros(3))
        else:
            # the two end particles are attached to one segment only and carry half the mass
            mass = 0.5 * m_particle if i in (0, n) else m_particle
            acc.append(se.g_earth + force[i] / mass)
    return acc


def build_model(se, fix_p1, fix_p2, v_ro):
    """ The tether as one CasADi expression graph. Returns the symbolic time `t`, the
        state vector `y` = (pos, vel) and its derivative `ydot`. """
    n = se.segments
    split = 3 * (n + 1)                        # start of the velocities in y
    t = ca.SX.sym('t')
    y = ca.SX.sym('y', 2 * split)
    pos = [y[3*i:3*i + 3] for i in range(n + 1)]
    vel = [y[split + 3*i:split + 3*i + 3] for i in range(n + 1)]
    acc = accelerations(t, pos, vel, se, fix_p1, fix_p2, v_ro)
    return t, y, ca.vertcat(*vel, *acc)


def find_steady_state(se, fix_p1, fix_p2, pos0):
    """ Find the steady-state tether shape for v_ro = 0 with scipy.optimize.least_squares:
        solve for the positions of the non-fixed particles such that their acceleration
        (with all velocities zero) is zero. `pos0` supplies the initial guess and the
        (unchanged) positions of the fixed endpoint(s). CasADi supplies the Jacobian.

        least_squares (trust-region, x_scale='jac') is used rather than root(method='hybr'):
        the huge spread between c_spring's magnitude and typical position values makes the
        problem badly scaled, and hybr's fixed internal scaling fails to converge on it,
        while least_squares' automatic Jacobian-based scaling handles it reliably. """
    n = se.segments
    split = 3 * (n + 1)
    free_idx = [i for i in range(n + 1) if not ((i == 0 and fix_p1) or (i == n and fix_p2))]
    free_dof = np.concatenate([np.arange(3 * i, 3 * i + 3) for i in free_idx])

    t, y, ydot = build_model(se, fix_p1, fix_p2, 0.0)
    acc_sym = ydot[split:]
    f_acc = ca.Function('acc', [y], [acc_sym])
    f_jac = ca.Function('acc_jac', [y], [ca.jacobian(acc_sym, y[:split])])

    # start from a sagging shape, not from the straight line, so that the solver converges
    # to the stable hanging equilibrium rather than to the unstable arch (see add_initial_sag)
    x0 = add_initial_sag(se, pos0)[free_idx].flatten()

    def state(x):
        pos = pos0.copy()
        pos[free_idx] = x.reshape(len(free_idx), 3)
        return np.concatenate([pos.flatten(), np.zeros(split)])

    def residual(x):
        return np.array(f_acc(state(x))).flatten()[free_dof]

    def jacobian(x):
        return np.array(f_jac(state(x)))[np.ix_(free_dof, free_dof)]

    # the axial spring stiffness is many orders of magnitude larger than the gravity/drag
    # forces that bend the tether out of a straight line, which makes this system badly
    # scaled and needs a generous max_nfev; x_scale='jac' keeps it tractable.
    sol = least_squares(residual, x0, jac=jacobian, method='trf', x_scale='jac',
                        xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=200_000)
    if not sol.success or np.max(np.abs(sol.fun)) > 1e-6:
        raise RuntimeError(f"Steady state solver failed: {sol.message}")
    pos = pos0.copy()
    pos[free_idx] = sol.x.reshape(len(free_idx), 3)
    return pos


def simulate(se, pos0, vel0, fix_p1, fix_p2):
    """ Simulate the tether model from the initial condition (pos0, vel0) over the
        duration se.duration with CVODES and CasADi's analytic sparse Jacobian, storing
        the result on a 0.02 s grid. Returns (t_sol, y, elapsed_time). """
    dt = 0.02
    t, y, ydot = build_model(se, fix_p1, fix_p2, se.v_ro)
    y0 = np.concatenate([pos0.flatten(), vel0.flatten()])
    t_sol = np.linspace(0.0, se.duration, round(se.duration / dt) + 1)
    # 'csparse' factorises the Jacobian CasADi derived from the model; its sparsity is
    # found by CasADi and is block-tridiagonal, so the dense solver would do most of its
    # work on structural zeros
    sim = ca.integrator('sim', 'cvodes', {'x': y, 't': t, 'ode': ydot}, 0.0, t_sol[1:],
                        {'abstol': 1.0e-6, 'reltol': 1.0e-6,
                         'linear_multistep_method': 'bdf',
                         'nonlinear_solver_iteration': 'newton',
                         'linear_solver': 'csparse'})
    start = time_module.time()
    xf = np.array(sim(x0=y0)['xf'])
    elapsed = time_module.time() - start
    return t_sol, np.column_stack([y0, xf]).T, elapsed


def plot2d(fig, pos, t, se, line, sc, txt):
    x, z = pos[:, 0], pos[:, 2]
    if line is None:
        line, = plt.plot(x, z, linewidth=1)
        sc = plt.scatter(x, z, s=15, color="red")
        txt = plt.annotate(f"t={t:.1f} s", xy=(se.l0 / 4.2, -7.0), fontsize=12)
        plt.show(block=False)
    else:
        line.set_xdata(x)
        line.set_ydata(z)
        sc.set_offsets(np.c_[x, z])
        txt.set_text(f"t={t:.1f} s")
        fig.canvas.draw()
    plt.pause(0.01)
    return line, sc, txt


def play(se, t_sol, y):
    """ Animate the solution, plotting the tether shape every 151 ms of simulated time,
        at half of real-time speed. If se.save is True, one PNG per frame is written to
        the video folder. """
    n = se.segments
    dt = 0.151
    plt.ioff()
    fig = plt.figure()
    plt.ylim(-1.2 * (se.l0 + se.v_ro * se.duration), 0.5)
    plt.xlim(-se.l0, se.l0)
    plt.grid(True, color="grey", linestyle="dotted")
    if se.save:
        os.makedirs("video", exist_ok=True)
    line, sc, txt = None, None, None
    j = 0
    start = time_module.time()
    for t in np.arange(0.0, se.duration, dt):
        idx = min(np.searchsorted(t_sol, t), t_sol.shape[0] - 1)
        pos = y[idx, :3 * (n + 1)].reshape(n + 1, 3)
        line, sc, txt = plot2d(fig, pos, t, se, line, sc, txt)
        if se.save:
            plt.savefig(f"video/img-{j:04d}.png")
        j += 1
        wait_until = start + 0.5 * t
        now = time_module.time()
        if wait_until > now:
            time_module.sleep(wait_until - now)
    if se.save:
        print("Run the script ./bin/export_gif to create the gif file!")
    if os.environ.get("TETHERS_BRIEF_PLOT") == "1":
        plt.pause(1)
        plt.close('all')
    else:
        plt.show(block=True)


def main(p1=np.array([0.0, 0.0, 0.0]), p2=None, fix_p1=True, fix_p2=False):
    """ Build, simulate and animate the tether model for the endpoints p1 and p2 (see
        calc_initial_state / find_steady_state for their meaning).
        Returns (t_sol, y, se). """
    se = Settings3()
    set_tether_diameter(se, se.d_tether)  # adapt spring and damping constants to tether diameter
    pos0, vel0, p1, p2 = calc_initial_state(se, p1, p2)
    v_ro = se.v_ro          # save the reel-out speed
    se.v_ro = 0.0           # v_ro must be zero to find the steady state
    try:
        # the steady-state shape is always found with BOTH endpoints held fixed at their
        # straight-line positions (matching Tether_08.jl's `model`, which always solves the
        # steady state with fix_p1=true, fix_p2=true), regardless of the fix_p1/fix_p2 the
        # caller wants for the actual simulation; those are applied only afterwards.
        pos_ss = find_steady_state(se, True, True, pos0)
    finally:
        se.v_ro = v_ro       # restore the reel-out speed, also if the solver failed
    t_sol, y, elapsed_time = simulate(se, pos_ss, vel0, fix_p1, fix_p2)

    # saving the z position and velocity of the last (free) particle for comparison
    # with the Julia implementation
    n = se.segments
    pos_z = y[:, 3 * n + 2]
    vel_z = y[:, 3 * (n + 1) + 3 * n + 2]
    os.makedirs("output", exist_ok=True)
    with open(os.path.join("output", "Tether_08_python.csv"), "w") as f:
        f.write("time,pos_z,vel_z\n")
        for t_i, pz_i, vz_i in zip(t_sol, pos_z, vel_z):
            f.write(f"{t_i},{pz_i},{vz_i}\n")

    if os.environ.get("TETHERS_PRECOMPILE") == "1":
        return t_sol, y, se
    play(se, t_sol, y)
    print(f"Elapsed time: {elapsed_time} s, speed: {round(se.duration / elapsed_time)} times real-time")
    return t_sol, y, se


if __name__ == '__main__':
    # run the simulation with a free (unfixed) second attachment point
    main(p2=np.array([-40.0, 0.0, -47.0]), fix_p2=False)
