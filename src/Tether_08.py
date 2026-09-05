# -*- coding: utf-8 -*-
"""
Tutorial example simulating a 3D mass-spring system with a nonlinear spring (1% stiffness
for l < l_0), n tether segments, tether drag and reel-in and reel-out.

New feature: instead of ModelingToolkit's symbolic steady-state machinery (used by
Tether_08.jl), the initial tether shape for a given pair of endpoints is found directly
with scipy.optimize.root, solving for the positions at which the acceleration of every
free particle is zero (with velocities zero and v_ro = 0). That steady-state shape is
then used as the initial condition for the time simulation.
"""
import os
import math
import time as time_module

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
import matplotlib.pyplot as plt

from dataclasses import dataclass, field


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
        Returns (pos0, vel0, p1, p2), each pos/vel a 3 x (segments+1) array. """
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
    pos0 = np.zeros((3, n + 1))
    vel0 = np.zeros((3, n + 1))
    delta = (p2 - p1) / n
    for i in range(n + 1):
        pos0[:, i] = p1 + i * delta
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
    chord = np.linalg.norm(pos0[:, n] - pos0[:, 0])
    if chord >= se.l0:
        return pos0                     # the tether is taut, the straight line is fine
    # parabolic bow: zero at both endpoints, maximum in the middle
    xi = np.linspace(0.0, 1.0, n + 1)
    shape = 4.0 * xi * (1.0 - xi)

    def polyline_length(sag):
        pos = pos0.copy()
        pos[2, :] -= sag * shape
        return np.sum(np.linalg.norm(pos[:, 1:] - pos[:, :-1], axis=0))

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
    pos[2, :] -= 0.5 * (lo + hi) * shape
    return pos


def calc_accelerations(t, pos, vel, se, fix_p1, fix_p2, v_ro):
    """ Acceleration of every tether particle, given the current positions `pos` and
        velocities `vel` (each a 3 x (segments+1) array), the reel-out speed `v_ro` and
        whether the two endpoints are held fixed. Mirrors the equations built by
        Tether_08.jl's `model` function. """
    n = se.segments
    l_spring = (se.l0 + v_ro * t) / n
    c_spring = se.c_spring / l_spring
    damping = se.damping / l_spring
    mass_per_meter = se.rho_tether * math.pi * (se.d_tether / 2000.0) ** 2
    m_tether_particle = mass_per_meter * l_spring

    segment = pos[:, 1:] - pos[:, :-1]                      # (3, n)
    length = np.linalg.norm(segment, axis=0)                # (n,)
    unit_vector = -segment / length                         # (3, n), points from i+1 towards i
    rel_vel = vel[:, 1:] - vel[:, :-1]                       # (3, n)
    spring_vel = -np.einsum('ij,ij->j', unit_vector, rel_vel)
    rcs = se.rel_compression_stiffness
    c_spr = c_spring / (1.0 + rcs) * (rcs + (length > l_spring))
    spring_force = (c_spr * (length - l_spring) + damping * spring_vel) * unit_vector

    v_apparent = se.v_wind_tether[:, None] - (vel[:, :-1] + vel[:, 1:]) / 2.0
    v_app_along = np.einsum('ij,ij->j', v_apparent, unit_vector)
    v_app_perp = v_apparent - v_app_along * unit_vector
    norm_v_app = np.linalg.norm(v_app_perp, axis=0)
    half_drag_force = (0.25 * se.rho * se.cd_tether * norm_v_app
                        * (length * se.d_tether / 1000.0)) * v_app_perp

    acc = np.zeros_like(pos)
    for i in range(n + 1):
        if i == n:
            total_force = spring_force[:, i - 1] + half_drag_force[:, i - 1]
            acc[:, i] = np.zeros(3) if fix_p2 else se.g_earth + total_force / (0.5 * m_tether_particle)
        elif i == 0:
            total_force = -spring_force[:, i] + half_drag_force[:, i]
            acc[:, i] = np.zeros(3) if fix_p1 else se.g_earth + total_force / (0.5 * m_tether_particle)
        else:
            total_force = (spring_force[:, i - 1] - spring_force[:, i]
                            + half_drag_force[:, i - 1] + half_drag_force[:, i])
            acc[:, i] = se.g_earth + total_force / m_tether_particle
    return acc


def find_steady_state(se, fix_p1, fix_p2, pos0):
    """ Find the steady-state tether shape for v_ro = 0 with scipy.optimize.least_squares:
        solve for the positions of the non-fixed particles such that their acceleration
        (with all velocities zero) is zero. `pos0` supplies the initial guess and the
        (unchanged) positions of the fixed endpoint(s).

        least_squares (trust-region, x_scale='jac') is used rather than root(method='hybr'):
        the huge spread between c_spring's magnitude and typical position values makes the
        problem badly scaled, and hybr's fixed internal scaling fails to converge on it,
        while least_squares' automatic Jacobian-based scaling handles it reliably. """
    n = se.segments
    free_idx = [i for i in range(n + 1) if not ((i == 0 and fix_p1) or (i == n and fix_p2))]
    # start from a sagging shape, not from the straight line, so that the solver converges
    # to the stable hanging equilibrium rather than to the unstable arch (see add_initial_sag)
    x0 = add_initial_sag(se, pos0)[:, free_idx].flatten()

    def residual(x):
        pos = pos0.copy()
        pos[:, free_idx] = x.reshape(3, len(free_idx))
        vel = np.zeros_like(pos)
        acc = calc_accelerations(0.0, pos, vel, se, fix_p1, fix_p2, 0.0)
        return acc[:, free_idx].flatten()

    # the axial spring stiffness is many orders of magnitude larger than the gravity/drag
    # forces that bend the tether out of a straight line, which makes this system badly
    # scaled and needs a generous max_nfev; x_scale='jac' keeps it tractable.
    sol = least_squares(residual, x0, method='trf', x_scale='jac',
                         xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=200_000)
    if not sol.success or np.max(np.abs(sol.fun)) > 1e-6:
        raise RuntimeError(f"Steady state solver failed: {sol.message}")
    pos_ss = pos0.copy()
    pos_ss[:, free_idx] = sol.x.reshape(3, len(free_idx))
    return pos_ss


def _rhs(t, state, se, fix_p1, fix_p2, v_ro):
    n = se.segments
    m = 3 * (n + 1)
    pos = state[:m].reshape(3, n + 1)
    vel = state[m:].reshape(3, n + 1)
    acc = calc_accelerations(t, pos, vel, se, fix_p1, fix_p2, v_ro)
    return np.concatenate([vel.flatten(), acc.flatten()])


def simulate(se, pos0, vel0, fix_p1, fix_p2):
    """ Simulate the tether model from the initial condition (pos0, vel0) over the
        duration se.duration with an implicit multistep solver (BDF, SciPy's counterpart
        to OrdinaryDiffEq's FBDF), storing the result on a 0.02 s grid.
        Returns (sol, elapsed_time). """
    dt = 0.02
    t_eval = np.arange(0.0, se.duration + 1e-9, dt)
    state0 = np.concatenate([pos0.flatten(), vel0.flatten()])
    start = time_module.time()
    sol = solve_ivp(_rhs, (0.0, se.duration), state0, method='BDF', t_eval=t_eval,
                     args=(se, fix_p1, fix_p2, se.v_ro), rtol=1e-6, atol=1e-6)
    elapsed = time_module.time() - start
    if not sol.success:
        raise RuntimeError(f"Simulation failed: {sol.message}")
    return sol, elapsed


def plot2d(fig, pos, t, se, line, sc, txt):
    x, z = pos[0, :], pos[2, :]
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


def play(se, sol):
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
    m = 3 * (n + 1)
    start = time_module.time()
    for t in np.arange(0.0, se.duration, dt):
        idx = min(np.searchsorted(sol.t, t), sol.t.shape[0] - 1)
        pos = sol.y[:m, idx].reshape(3, n + 1)
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
        calc_initial_state / find_steady_state for their meaning). Returns (sol, se). """
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
    sol, elapsed_time = simulate(se, pos_ss, vel0, fix_p1, fix_p2)

    # saving the z position and velocity of the last (free) particle for comparison
    # with the Julia implementation
    n = se.segments
    pos_z = sol.y[3 * n + 2, :]
    vel_z = sol.y[3 * (n + 1) + 3 * n + 2, :]
    os.makedirs("output", exist_ok=True)
    with open(os.path.join("output", "Tether_08_python.csv"), "w") as f:
        f.write("time,pos_z,vel_z\n")
        for t_i, pz_i, vz_i in zip(sol.t, pos_z, vel_z):
            f.write(f"{t_i},{pz_i},{vz_i}\n")

    if os.environ.get("TETHERS_PRECOMPILE") == "1":
        return sol, se
    play(se, sol)
    print(f"Elapsed time: {elapsed_time} s, speed: {round(se.duration / elapsed_time)} times real-time")
    return sol, se


if __name__ == '__main__':
    # run the simulation with a free (unfixed) second attachment point
    main(p2=np.array([-40.0, 0.0, -47.0]), fix_p2=False)
