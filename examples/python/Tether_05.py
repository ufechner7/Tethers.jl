# -*- coding: utf-8 -*-
"""
Tutorial example simulating a 3D mass-spring system with a nonlinear spring (no spring
force while a segment is loose) and five tether segments.

The model is written as a CasADi expression graph and integrated with SUNDIALS' CVODES,
which CasADi supplies with the exact Jacobian of the model.

State vector y = pos0, pos1, vel1, ..., posN, velN (each a 3D vector). pos0 is the fixed
attachment point, so its derivative is zero and it stays where it starts.
"""
import math
import os
import time

import numpy as np
import matplotlib.pyplot as plt
import casadi as ca

G_EARTH  = np.array([0.0, 0.0, -9.81]) # gravitational acceleration
C_SPRING = 50.0                        # spring constant
DAMPING  =  0.5                        # damping [Ns/m]
L0      =  5.0                         # initial segment length     [m]
V0       =  2.0                        # initial velocity of the lowest mass [m/s]
ALPHA0   = math.pi/10                  # initial tether angle     [rad]
SEGMENTS = 5
MASS     = 0.5                         # mass per tether particle  [kg]
DURATION = 10                          # duration of the simulation [s]
NONLINEAR = True                       # if True, a loose segment exerts no spring force


def calc_spring_constant(norm):
    """ Spring constant of one segment. A loose segment (norm <= L0) cannot push,
        therefore its spring constant is zero if NONLINEAR is True. """
    if not NONLINEAR:
        return C_SPRING
    return ca.if_else(norm > L0, C_SPRING, 0.0)


def calc_spring_force(pos1, pos2, vel1, vel2):
    """ Spring and damping force of the segment between the masses at pos1 and pos2.
        The result points from pos1 towards pos2. Spring and damper act in parallel
        along the segment, therefore the damping force uses the component of the
        relative velocity along the segment and not the full 3D vector. """
    segment     = pos2 - pos1
    norm        = ca.norm_2(segment)
    unit_vector = segment / norm
    spring_vel  = ca.dot(vel2 - vel1, unit_vector) # rate of change of the segment length
    return (calc_spring_constant(norm) * (norm - L0) + DAMPING * spring_vel) * unit_vector


def calc_particle_mass(i):
    """ Mass of particle i. The last particle is attached to one segment only and
        therefore carries half of the mass of an inner particle. """
    if i == SEGMENTS:
        return 0.5 * MASS
    return MASS


def pos_block(k):
    """ Index, in 3-vector blocks of the state vector, of the position of mass k. """
    return 0 if k == 0 else 2*k - 1


def vel_block(k):
    """ Index, in 3-vector blocks of the state vector, of the velocity of mass k (k >= 1). """
    return 2*k


def build_model():
    """ The tether as one CasADi expression graph. Returns the state vector `y`, its
        derivative `ydot` and the initial state `y0`. """
    n = SEGMENTS
    y = ca.SX.sym('y', 3 * (2*n + 1))
    blocks = [y[3*i:3*i + 3] for i in range(2*n + 1)]
    pos = [blocks[pos_block(k)] for k in range(n + 1)]
    vel = [ca.SX.zeros(3)] + [blocks[vel_block(k)] for k in range(1, n + 1)]

    force = [calc_spring_force(pos[k], pos[k+1], vel[k], vel[k+1]) for k in range(n)]
    ydot = [ca.SX.zeros(3)] * (2*n + 1)
    for k in range(1, n + 1):
        # no segment below the last particle
        below = force[k] if k < n else ca.SX.zeros(3)
        ydot[pos_block(k)] = vel[k]
        ydot[vel_block(k)] = G_EARTH + (below - force[k-1]) / calc_particle_mass(k)

    y0 = np.zeros(3 * (2*n + 1))
    for k in range(n + 1):
        l, v = -k * L0, k * V0 / n
        y0[3*pos_block(k):3*pos_block(k) + 3] = [math.sin(ALPHA0) * l, 0.0, math.cos(ALPHA0) * l]
        if k > 0:
            y0[3*vel_block(k):3*vel_block(k) + 3] = [math.sin(ALPHA0) * v, 0.0, math.cos(ALPHA0) * v]
    return y, ca.vertcat(*ydot), y0


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
    dt = 0.15
    plt.ioff()
    fig = plt.figure()
    plt.ylim(-SEGMENTS*L0-10, 0.5)
    plt.xlim(-SEGMENTS*L0/2, SEGMENTS*L0/2)
    plt.grid(True, color="grey", linestyle="dotted")
    line, sc, txt = None, None, None
    for t in np.linspace(0, duration, num=round(duration/dt)+1):
        line, sc, txt = plot2d(fig, t_sol, y, t, SEGMENTS, line, sc, txt)
        time.sleep(dt/2)
    if os.environ.get("TETHERS_BRIEF_PLOT") == "1":
        # show briefly and close automatically, e.g. when running the tests
        plt.pause(1)
        plt.close('all')
    else:
        plt.show()

def run_example():
    y, ydot, y0 = build_model()
    t_sol = np.linspace(0.0, DURATION, round(DURATION*50) + 1)   # 50 points per second
    sim = ca.integrator('sim', 'cvodes', {'x': y, 'ode': ydot}, 0.0, t_sol[1:],
                        {'abstol': 1.0e-6, 'reltol': 1.0e-6})
    y_sol = np.column_stack([y0, np.array(sim(x0=y0)['xf'])]).T

    # extract the z position and velocity of the lowest mass (mass SEGMENTS)
    pos_z_ix = 5 + (SEGMENTS - 1) * 6
    vel_z_ix = pos_z_ix + 3
    pos_z = y_sol[:, pos_z_ix]
    vel_z = y_sol[:, vel_z_ix]

    # saving the result for comparison with the Julia implementation
    os.makedirs("output", exist_ok=True)
    with open(os.path.join("output", "Tether_05_python.csv"), "w") as f:
        f.write("time,pos_z,vel_z\n")
        for t_i, pz_i, vz_i in zip(t_sol, pos_z, vel_z):
            f.write(f"{t_i},{pz_i},{vz_i}\n")

    play(DURATION, t_sol, y_sol)


if __name__ == '__main__':
    run_example()
