# -*- coding: utf-8 -*-
"""
Tutorial example simulating a 3D mass-spring system with a nonlinear spring (1%
stiffness for l < l_0), five tether segments, tether drag and reel-out.

The model is written once as a CasADi expression graph. CasADi differentiates it to get
the exact Jacobian and finds its sparsity, and SUNDIALS' CVODES integrates it with a BDF
formula and a sparse Newton solve - the same fixed-leading-coefficient BDF that the FBDF
solver of Tether_07.jl uses. This is the Python counterpart of Tether_07.jl, extending
Tether_06.py with an aerodynamic drag force per segment.

State vector y = pos0, pos1, vel1, ..., posN, velN (each a 3D vector). pos0 is the fixed
attachment point, so its derivative is zero and it stays where it starts.
"""
import math
import os
import time

import numpy as np
import matplotlib.pyplot as plt
import casadi as ca

G_EARTH       = np.array([0.0, 0.0, -9.81])   # gravitational acceleration        [m/s²]
V_WIND_TETHER = np.array([2.0, 0.0, 0.0])     # wind velocity acting on the tether  [m/s]
RHO           = 1.225                          # air density                     [kg/m³]
CD_TETHER     = 0.958                          # drag coefficient of the tether       [-]
C_SPRING      =  614600.0                      # spring constant
DAMPING       =  473                           # unit damping constant [Ns], must match Tether_07.jl
REL_COMPRESSION_STIFFNESS = 0.01               # relative compression stiffness       [-]
L0            =  50.0                          # initial segment length     [m]
ALPHA0        = math.pi/10                     # initial tether angle     [rad]
SEGMENTS      = 5
DURATION      = 10                             # duration of the simulation [s]

V_RO = 2.0                             # reel-out speed                  [m/s]
D_TETHER = 4.0                         # tether diameter                  [mm]
RHO_TETHER = 724.0                     # density of Dyneema            [kg/m³]

mass_per_meter = RHO_TETHER * math.pi * (D_TETHER/2000.0)**2


def calc_c_spring(norm, l_seg, c_spring0):
    """ Nonlinear spring stiffness of a segment: full stiffness while taut
        (norm > l_seg), only REL_COMPRESSION_STIFFNESS of it while slack. Unlike
        Tether_06.py, the switch is a hard step, matching Tether_07.jl. `if_else` is
        differentiated on the active branch, so the Jacobian is the one of the branch the
        segment is currently on, and the step itself contributes nothing. """
    taut = ca.if_else(norm > l_seg, 1.0, 0.0)
    return c_spring0 / (1.0 + REL_COMPRESSION_STIFFNESS) * (REL_COMPRESSION_STIFFNESS + taut)


def calc_spring_force(pos1, pos2, vel1, vel2, l_seg, c_spring0, damping):
    """ Spring and damping force of the segment between the masses at pos1 and pos2.
        The result points from pos1 towards pos2. Spring and damper act in parallel
        along the segment, therefore the damping force uses the component of the
        relative velocity along the segment and not the full 3D vector. """
    segment     = pos2 - pos1
    norm        = ca.norm_2(segment)
    unit_vector = segment / norm
    c_spring    = calc_c_spring(norm, l_seg, c_spring0)
    spring_vel  = ca.dot(vel2 - vel1, unit_vector)   # rate of change of the segment length
    return (c_spring * (norm - l_seg) + damping * spring_vel) * unit_vector


def calc_drag_force(pos1, pos2, vel1, vel2):
    """ Half of the aerodynamic drag force of the segment between pos1 and pos2,
        matching half_drag_force in Tether_07.jl. Unlike the spring force, this is
        not an action/reaction pair: the same vector is added to both end masses,
        since it is one half of the drag of the segment they share. """
    segment     = pos2 - pos1
    norm        = ca.norm_2(segment)
    unit_vector = segment / norm
    v_apparent  = V_WIND_TETHER - (vel1 + vel2) / 2.0
    v_app_perp  = v_apparent - ca.dot(v_apparent, unit_vector) * unit_vector
    norm_v_app  = ca.norm_2(v_app_perp)
    return 0.25 * RHO * CD_TETHER * norm_v_app * (norm * D_TETHER / 1000.0) * v_app_perp


def calc_particle_mass(i, m_tether_particle):
    """ Mass of particle i. The last particle is attached to one segment only and
        therefore carries half of the mass of an inner particle. """
    if i == SEGMENTS:
        return 0.5 * m_tether_particle
    return m_tether_particle


def pos_block(k):
    """ Index, in 3-vector blocks of the state vector, of the position of mass k. """
    return 0 if k == 0 else 2*k - 1


def vel_block(k):
    """ Index, in 3-vector blocks of the state vector, of the velocity of mass k (k >= 1). """
    return 2*k


def build_model():
    """ The tether as one CasADi expression graph. Returns the symbolic time `t`, the
        state vector `y`, its derivative `ydot` and the initial state `y0`. """
    n = SEGMENTS
    t = ca.SX.sym('t')
    y = ca.SX.sym('y', 3 * (2*n + 1))
    blocks = [y[3*i:3*i + 3] for i in range(2*n + 1)]
    pos = [blocks[pos_block(k)] for k in range(n + 1)]
    vel = [ca.SX.zeros(3)] + [blocks[vel_block(k)] for k in range(1, n + 1)]

    l_seg = (L0 + V_RO*t) / n          # unstretched length of one segment
    c_spring0 = C_SPRING / l_seg
    damping = DAMPING / l_seg
    m_tether_particle = mass_per_meter * l_seg

    force = [calc_spring_force(pos[k], pos[k+1], vel[k], vel[k+1],
                               l_seg, c_spring0, damping) for k in range(n)]
    drag  = [calc_drag_force(pos[k], pos[k+1], vel[k], vel[k+1]) for k in range(n)]

    ydot = [ca.SX.zeros(3)] * (2*n + 1)
    ydot[0] = ca.SX.zeros(3)           # mass 0 is fixed at the attachment point
    for k in range(1, n + 1):
        zero = ca.SX.zeros(3)
        # the spring force is an action/reaction pair, so the two segments of a mass
        # cancel each other, while both segments push it downwind with their drag half
        spring_forces = force[k-1] - (force[k] if k < n else zero)
        drag_forces   = drag[k-1] + (drag[k] if k < n else zero)
        mass = calc_particle_mass(k, m_tether_particle)
        ydot[pos_block(k)] = vel[k]
        ydot[vel_block(k)] = G_EARTH - (spring_forces - drag_forces) / mass

    y0 = np.zeros(3 * (2*n + 1))
    for k in range(n + 1):
        l, v = -k * L0 / n, -k * V_RO / n
        y0[3*pos_block(k):3*pos_block(k) + 3] = [math.sin(ALPHA0) * l, 0.0, math.cos(ALPHA0) * l]
        if k > 0:
            y0[3*vel_block(k):3*vel_block(k) + 3] = [math.sin(ALPHA0) * v, 0.0, math.cos(ALPHA0) * v]
    return t, y, ca.vertcat(*ydot), y0


def plot2d(fig, t_sol, y, reltime, segments, line, sc, txt):
    index = min(np.searchsorted(t_sol, reltime), len(t_sol) - 1)
    x, z = np.zeros(segments+1), np.zeros(segments+1)
    for i in range(segments):
        x[i+1] = y[index, 3+6*i]
        z[i+1] = y[index, 5+6*i]
    if line is None:
        z_max = np.max(z)
        line, = plt.plot(x, z, linewidth=1)
        sc  = plt.scatter(x, z, s=15, color="red")
        plt.pause(0.01)
        txt = plt.annotate("t="+str(round(reltime,1))+" s",
                           xy=(segments*L0/4.2, z_max-3.0*segments/5), fontsize = 12)
        plt.show(block=False)
    else:
        line.set_xdata(x)
        line.set_ydata(z)
        sc.set_offsets(np.c_[x, z])
        txt.set_text("t="+str(round(reltime,1))+" s")
        fig.canvas.draw()
        plt.pause(0.01)
        plt.show(block=False)
    return line, sc, txt

def play(duration, t_sol, y):
    dt = 0.151
    plt.ioff()
    fig = plt.figure()
    plt.ylim(-1.2*(L0+V_RO*duration), 0.5)
    plt.xlim(-L0/2, L0/2)
    plt.grid(True, color="grey", linestyle="dotted")
    plt.tight_layout(rect=(0, 0, 0.98, 0.98))
    line, sc, txt = None, None, None
    for t in np.linspace(0, duration, num=round(duration/dt)):
        line, sc, txt = plot2d(fig, t_sol, y, t, SEGMENTS, line, sc, txt)
        time.sleep(dt/2)
    if os.environ.get("TETHERS_BRIEF_PLOT") == "1":
        # show briefly and close automatically, e.g. when running the tests
        plt.pause(1)
        plt.close('all')
    else:
        plt.show(block=True)

def run_example():
    t, y, ydot, y0 = build_model()
    t_sol = np.linspace(0.0, DURATION, round(DURATION*50) + 1)   # 50 points per second
    # 'csparse' factorises the Jacobian CasADi derived from the model; its sparsity is
    # found by CasADi and is block-tridiagonal, so the dense solver would do most of its
    # work on structural zeros
    sim = ca.integrator('sim', 'cvodes', {'x': y, 't': t, 'ode': ydot}, 0.0, t_sol[1:],
                        {'abstol': 1.0e-6, 'reltol': 1.0e-6, 'max_multistep_order': 3,
                         'linear_multistep_method': 'bdf',
                         'nonlinear_solver_iteration': 'newton',
                         'linear_solver': 'csparse'})
    start = time.perf_counter()
    xf = np.array(sim(x0=y0)['xf'])
    elapsed_time = time.perf_counter() - start
    print(f"Elapsed time: {elapsed_time} s, speed: {round(DURATION/elapsed_time)} times real-time")
    y_sol = np.column_stack([y0, xf]).T

    # extract the z position and velocity of the lowest mass (mass SEGMENTS)
    pos_z_ix = 5 + (SEGMENTS - 1) * 6
    vel_z_ix = pos_z_ix + 3
    pos_z = y_sol[:, pos_z_ix]
    vel_z = y_sol[:, vel_z_ix]

    # saving the result for comparison with the Julia implementation
    os.makedirs("output", exist_ok=True)
    with open(os.path.join("output", "Tether_07_python.csv"), "w") as f:
        f.write("time,pos_z,vel_z\n")
        for t_i, pz_i, vz_i in zip(t_sol, pos_z, vel_z):
            f.write(f"{t_i},{pz_i},{vz_i}\n")

    play(DURATION, t_sol, y_sol)
    return


if __name__ == '__main__':
    run_example()
