# -*- coding: utf-8 -*-
"""
Tutorial example simulating a 3D mass-spring system with a nonlinear spring (1% stiffness
for l < l_0), n tether segments, tether drag and reel-in and reel-out.

New feature: instead of ModelingToolkit's symbolic steady-state machinery (used by
Tether_08.jl), the initial tether shape for a given pair of endpoints is found directly
with scipy.optimize.least_squares, solving for the positions at which the acceleration of
every free particle is zero (with velocities zero and v_ro = 0). That steady-state shape is
then used as the initial condition for the time simulation.

The time simulation uses Assimulo's IDA, the same solver as Tether_06.py and Tether_07.py,
with an analytic Jacobian. Both the steady-state solver and IDA are handed the analytic
derivatives of the particle accelerations w.r.t. the positions and velocities, which is
what makes this stiff system tractable: the spring stiffness is many orders of magnitude
larger than the gravity and drag forces, so finite differences are both slow (a residual
evaluation per state) and inaccurate here.

State vector y  = pos[0..n], vel[0..n]   (each particle contributes a 3D-vector)
Derivative   yd = vel[0..n], acc[0..n]
Residual    res = (yd.pos - vel), (yd.vel - acc)
"""
import os
import math
import time as time_module

import numpy as np
from scipy.optimize import least_squares
import matplotlib.pyplot as plt

from dataclasses import dataclass, field

from assimulo.solvers.sundials import IDA     # Imports the solver IDA from Assimulo
from assimulo.problem import Implicit_Problem # Imports the problem formulation from Assimulo


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


def calc_segment_props(se, t, v_ro):
    """ Unstretched segment length, segment spring constant, segment damping constant and
        the mass of an interior tether particle at the time t, for the reel-out speed
        `v_ro`. All of them change while the tether is reeled in or out. """
    n = se.segments
    l_spring = (se.l0 + v_ro * t) / n
    c_spring = se.c_spring / l_spring
    damping = se.damping / l_spring
    mass_per_meter = se.rho_tether * math.pi * (se.d_tether / 2000.0) ** 2
    return l_spring, c_spring, damping, mass_per_meter * l_spring


def calc_masses(se, m_tether_particle):
    """ Mass of each particle. The two end particles are attached to one segment only and
        therefore carry half of the mass of an interior particle. """
    mass = np.full(se.segments + 1, m_tether_particle)
    mass[0] = mass[se.segments] = 0.5 * m_tether_particle
    return mass


def calc_segment_forces(se, pos, vel, l_spring, c_spring, damping):
    """ For every segment j (between the particles j and j+1) the tension force `fs[j]`,
        pulling particle j towards particle j+1, and `drag[j]`, one half of the
        aerodynamic drag force of the segment, which is applied to both of its ends.

        The spring is nonlinear: a taut segment (length > l_spring) uses the full
        stiffness, a slack one only se.rel_compression_stiffness of it. Spring and damper
        act in parallel along the segment, therefore the damping force uses the component
        of the relative velocity along the segment and not the full 3D vector. """
    segment = pos[1:] - pos[:-1]                              # (n, 3)
    length = np.linalg.norm(segment, axis=1)                  # (n,)
    e = segment / length[:, None]                             # unit vector, j towards j+1
    rel_vel = vel[1:] - vel[:-1]
    spring_vel = np.einsum('ij,ij->i', e, rel_vel)            # rate of change of the length
    rcs = se.rel_compression_stiffness
    c_spr = c_spring / (1.0 + rcs) * (rcs + (length > l_spring))
    fs = (c_spr * (length - l_spring) + damping * spring_vel)[:, None] * e

    v_app = se.v_wind_tether - (vel[:-1] + vel[1:]) / 2.0     # apparent wind velocity
    v_app_perp = v_app - np.einsum('ij,ij->i', v_app, e)[:, None] * e
    norm_v_app = np.linalg.norm(v_app_perp, axis=1)
    drag = (0.25 * se.rho * se.cd_tether * norm_v_app
            * length * se.d_tether / 1000.0)[:, None] * v_app_perp
    return fs, drag


def calc_accelerations(t, pos, vel, se, fix_p1, fix_p2, v_ro):
    """ Acceleration of every tether particle, given the current positions `pos` and
        velocities `vel` (each a (segments+1) x 3 array), the reel-out speed `v_ro` and
        whether the two endpoints are held fixed. Mirrors the equations built by
        Tether_08.jl's `model` function. """
    n = se.segments
    l_spring, c_spring, damping, m_tether_particle = calc_segment_props(se, t, v_ro)
    fs, drag = calc_segment_forces(se, pos, vel, l_spring, c_spring, damping)

    force = np.zeros_like(pos)
    force[:-1] += fs + drag        # segment j pulls particle j towards particle j+1 ...
    force[1:] += -fs + drag        # ... and particle j+1 towards particle j
    acc = se.g_earth + force / calc_masses(se, m_tether_particle)[:, None]
    if fix_p1:
        acc[0] = 0.0
    if fix_p2:
        acc[n] = 0.0
    return acc


def calc_segment_force_jac(se, pa, pb, va, vb, l_spring, c_spring, damping):
    """ Tension force `fs` and half drag force `drag` of the segment between the particles
        at pa and pb (as in calc_segment_forces), plus their analytic Jacobians, each a
        3x3 matrix dF_i/dx_j:

        - `dfs_dpb` = d(fs)/d(pb) = -d(fs)/d(pa)
        - `dfs_dvb` = d(fs)/d(vb) = -d(fs)/d(va)
        - `dd_dpa`, `dd_dpb` = d(drag)/d(pa), d(drag)/d(pb)
        - `dd_dv`   = d(drag)/d(va) = d(drag)/d(vb); the drag depends on the mean of the
          two velocities only, therefore it is the same matrix for both.

        The stiffness jumps at length == l_spring (the nonlinear spring is a hard step,
        as in Tether_08.jl); that jump is ignored here, so the Jacobian is the one of the
        currently active branch. """
    eye = np.eye(3)
    s = pb - pa
    length = np.linalg.norm(s)
    e = s / length
    p_mat = (eye - np.outer(e, e)) / length     # d(e)/d(pb) = -d(e)/d(pa), symmetric

    rel_vel = vb - va
    spring_vel = e @ rel_vel
    rcs = se.rel_compression_stiffness
    c_spr = c_spring / (1.0 + rcs) * (rcs + (1.0 if length > l_spring else 0.0))
    f_mag = c_spr * (length - l_spring) + damping * spring_vel
    fs = f_mag * e
    # d(f_mag)/d(pb), a row vector: the length grows with e, the damping term with the
    # rotation of the segment (p_mat is symmetric, so rel_vel^T p_mat = p_mat @ rel_vel)
    d_fmag_dpb = c_spr * e + damping * (p_mat @ rel_vel)
    dfs_dpb = np.outer(e, d_fmag_dpb) + f_mag * p_mat
    dfs_dvb = damping * np.outer(e, e)

    v_app = se.v_wind_tether - (va + vb) / 2.0
    v_app_along = v_app @ e
    v_app_perp = v_app - v_app_along * e
    norm_v_app = np.linalg.norm(v_app_perp)
    k_drag = 0.25 * se.rho * se.cd_tether * se.d_tether / 1000.0
    drag = k_drag * length * norm_v_app * v_app_perp
    if norm_v_app > 0.0:
        # d(norm_v_app * v_app_perp)/d(v_app_perp)
        g_mat = norm_v_app * eye + np.outer(v_app_perp, v_app_perp) / norm_v_app
        # -d(v_app_perp)/d(pb) = d(v_app_perp)/d(pa)
        a_mat = np.outer(e, p_mat @ v_app) + v_app_along * p_mat
        # the drag grows with the segment length (d(length)/d(pb) = e) and changes with
        # the direction of the segment, which turns v_app_perp
        dd_dpb = k_drag * (np.outer(norm_v_app * v_app_perp, e) - length * (g_mat @ a_mat))
        dd_dpa = k_drag * (-np.outer(norm_v_app * v_app_perp, e) + length * (g_mat @ a_mat))
        # d(v_app_perp)/d(va) = d(v_app_perp)/d(vb) = -(I - e e^T)/2
        dd_dv = -0.5 * k_drag * length * (g_mat @ (eye - np.outer(e, e)))
    else:
        # d(norm_v_app * v_app_perp) vanishes together with v_app_perp
        dd_dpa = dd_dpb = dd_dv = np.zeros((3, 3))
    return fs, drag, dfs_dpb, dfs_dvb, dd_dpa, dd_dpb, dd_dv


def calc_acc_jac(t, pos, vel, se, fix_p1, fix_p2, v_ro):
    """ Analytic Jacobian of the accelerations returned by calc_accelerations, as the pair
        (d(acc)/d(pos), d(acc)/d(vel)), each of shape (3*(n+1), 3*(n+1)) with the particles
        in the same order as in the flattened `pos` and `vel` arrays. """
    n = se.segments
    l_spring, c_spring, damping, m_tether_particle = calc_segment_props(se, t, v_ro)
    mass = calc_masses(se, m_tether_particle)
    size = 3 * (n + 1)
    dacc_dpos = np.zeros((size, size))
    dacc_dvel = np.zeros((size, size))
    fixed = [(i == 0 and fix_p1) or (i == n and fix_p2) for i in range(n + 1)]

    def blk(i):
        return slice(3 * i, 3 * i + 3)

    for j in range(n):
        _, _, dfs_dpb, dfs_dvb, dd_dpa, dd_dpb, dd_dv = calc_segment_force_jac(
            se, pos[j], pos[j + 1], vel[j], vel[j + 1], l_spring, c_spring, damping)
        # the segment pushes/pulls particle j with (fs + drag) and particle j+1 with
        # (-fs + drag), so both get the same drag and opposite tension contributions
        for i, sign in ((j, 1.0), (j + 1, -1.0)):
            if fixed[i]:
                continue
            inv_m = 1.0 / mass[i]
            dacc_dpos[blk(i), blk(j + 1)] += inv_m * (sign * dfs_dpb + dd_dpb)
            dacc_dpos[blk(i), blk(j)] += inv_m * (-sign * dfs_dpb + dd_dpa)
            dacc_dvel[blk(i), blk(j + 1)] += inv_m * (sign * dfs_dvb + dd_dv)
            dacc_dvel[blk(i), blk(j)] += inv_m * (-sign * dfs_dvb + dd_dv)
    return dacc_dpos, dacc_dvel


def find_steady_state(se, fix_p1, fix_p2, pos0):
    """ Find the steady-state tether shape for v_ro = 0 with scipy.optimize.least_squares:
        solve for the positions of the non-fixed particles such that their acceleration
        (with all velocities zero) is zero. `pos0` supplies the initial guess and the
        (unchanged) positions of the fixed endpoint(s). The analytic Jacobian of
        calc_acc_jac is passed to the solver.

        least_squares (trust-region, x_scale='jac') is used rather than root(method='hybr'):
        the huge spread between c_spring's magnitude and typical position values makes the
        problem badly scaled, and hybr's fixed internal scaling fails to converge on it,
        while least_squares' automatic Jacobian-based scaling handles it reliably. """
    n = se.segments
    free_idx = [i for i in range(n + 1) if not ((i == 0 and fix_p1) or (i == n and fix_p2))]
    free_dof = np.concatenate([np.arange(3 * i, 3 * i + 3) for i in free_idx])
    # start from a sagging shape, not from the straight line, so that the solver converges
    # to the stable hanging equilibrium rather than to the unstable arch (see add_initial_sag)
    x0 = add_initial_sag(se, pos0)[free_idx].flatten()

    def positions(x):
        pos = pos0.copy()
        pos[free_idx] = x.reshape(len(free_idx), 3)
        return pos

    def residual(x):
        pos = positions(x)
        acc = calc_accelerations(0.0, pos, np.zeros_like(pos), se, fix_p1, fix_p2, 0.0)
        return acc[free_idx].flatten()

    def jacobian(x):
        pos = positions(x)
        dacc_dpos, _ = calc_acc_jac(0.0, pos, np.zeros_like(pos), se, fix_p1, fix_p2, 0.0)
        return dacc_dpos[np.ix_(free_dof, free_dof)]

    # the axial spring stiffness is many orders of magnitude larger than the gravity/drag
    # forces that bend the tether out of a straight line, which makes this system badly
    # scaled and needs a generous max_nfev; x_scale='jac' keeps it tractable.
    sol = least_squares(residual, x0, jac=jacobian, method='trf', x_scale='jac',
                        xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=200_000)
    if not sol.success or np.max(np.abs(sol.fun)) > 1e-6:
        raise RuntimeError(f"Steady state solver failed: {sol.message}")
    return positions(sol.x)


# Extend Assimulos problem definition
class ExtendedProblem(Implicit_Problem):
    """ The tether model as an implicit DAE res(t, y, yd) = 0 for Assimulo's IDA. All
        states are differential here: the fixed endpoints are not enforced with an
        algebraic constraint, their acceleration is simply set to zero (they start at
        rest, therefore they stay where they are). """

    def __init__(self, se, pos0, vel0, fix_p1, fix_p2):
        self.se, self.fix_p1, self.fix_p2 = se, fix_p1, fix_p2
        self.split = 3 * (se.segments + 1)          # start of the velocities in y
        acc0 = calc_accelerations(0.0, pos0, vel0, se, fix_p1, fix_p2, se.v_ro)
        Implicit_Problem.__init__(self,
                                  y0=np.concatenate([pos0.flatten(), vel0.flatten()]),
                                  yd0=np.concatenate([vel0.flatten(), acc0.flatten()]),
                                  t0=0.0)
        self.name = 'Tether with drag and reel-out'

    def _unpack(self, y):
        return y[:self.split].reshape(-1, 3), y[self.split:].reshape(-1, 3)

    def res(self, t, y, yd):
        pos, vel = self._unpack(y)
        posd, veld = self._unpack(yd)
        acc = calc_accelerations(t, pos, vel, self.se, self.fix_p1, self.fix_p2, self.se.v_ro)
        return np.concatenate([(posd - vel).flatten(), (veld - acc).flatten()])

    def jac(self, c, t, y, yd):
        """ Analytic Jacobian J = d(res)/dy + c * d(res)/dyd, replacing the
            finite-difference approximation IDA would otherwise use. """
        pos, vel = self._unpack(y)
        m = self.split
        eye = np.eye(m)
        dacc_dpos, dacc_dvel = calc_acc_jac(t, pos, vel, self.se, self.fix_p1,
                                            self.fix_p2, self.se.v_ro)
        j_mat = np.zeros((2 * m, 2 * m))
        j_mat[:m, :m] = c * eye              # d(posd - vel)/d(pos), via yd
        j_mat[:m, m:] = -eye                 # d(posd - vel)/d(vel)
        j_mat[m:, :m] = -dacc_dpos           # d(veld - acc)/d(pos)
        j_mat[m:, m:] = c * eye - dacc_dvel  # d(veld - acc)/d(vel), plus the yd part
        return j_mat


def simulate(se, pos0, vel0, fix_p1, fix_p2):
    """ Simulate the tether model from the initial condition (pos0, vel0) over the
        duration se.duration with Assimulo's IDA and an analytic Jacobian, storing the
        result on a 0.02 s grid. Returns (t_sol, y, elapsed_time). """
    dt = 0.02
    model = ExtendedProblem(se, pos0, vel0, fix_p1, fix_p2)
    sim = IDA(model)  # Create the solver, using the default dense direct linear solver
    sim.verbosity = 0
    sim.atol = 1.0e-6
    sim.rtol = 1.0e-6
    sim.algvar = np.ones(2 * model.split)   # all states are differential
    sim.usejac = True

    start = time_module.time()
    t_raw, y_raw, _ = sim.simulate(se.duration, round(se.duration / dt))
    elapsed = time_module.time() - start

    # IDA's chosen output points can silently land off the shared 0.02 s grid. Resample the
    # full state onto the exact grid Tether_08.jl uses (ts=0:dt:duration), so the two CSVs
    # are directly comparable and play() gets a consistent (t_sol, y) pair.
    t_sol = np.linspace(0.0, se.duration, round(se.duration / dt) + 1)
    y = np.column_stack([np.interp(t_sol, t_raw, y_raw[:, j]) for j in range(y_raw.shape[1])])
    return t_sol, y, elapsed


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
