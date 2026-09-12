# -*- coding: utf-8 -*-
"""
Tutorial example simulating three falling masses, connected with two springs with damping.

The model is written as a CasADi expression graph and integrated with SUNDIALS' CVODES,
which CasADi supplies with the exact Jacobian of the model.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import casadi as ca

G_EARTH  = np.array([0.0, 0.0, -9.81]) # gravitational acceleration
C_SPRING = 50.0                        # spring constant [N/m]
DAMPING  =  0.5                        # damping [Ns/m]
L_0      = 10.0                        # initial segment length [m]
V0       = 4.0                         # initial velocity of the lowest mass [m/s]
SEGMENTS = 2                           # number of tether segments, must match Tether_04.jl
MASS     = 1.0                         # mass per tether of initial segment length [kg]


def pos_block(k):
    """ Index, in 3-vector blocks of the state vector, of the position of mass k. """
    return 0 if k == 0 else 2*k - 1


def vel_block(k):
    """ Index, in 3-vector blocks of the state vector, of the velocity of mass k (k >= 1). """
    return 2*k


def calc_spring_force(pos1, pos2, vel1, vel2):
    """ Spring and damping force of the segment between the masses at pos1 and pos2,
        pointing from pos1 towards pos2. Unlike the later examples the damper here acts
        on the full relative velocity vector, not only on its component along the
        segment. A slack segment (norm <= L_0) exerts no spring force. """
    segment = pos2 - pos1
    norm = ca.norm_2(segment)
    c_spring = ca.if_else(norm > L_0, C_SPRING, 0.0)
    return c_spring * (norm - L_0) * segment / norm + DAMPING * (vel2 - vel1)


def build_model():
    """ The chain of masses as one CasADi expression graph. Mass 0 is fixed, so its
        derivative is zero. Returns the state vector `y`, its derivative `ydot` and the
        initial state `y0`. """
    n = SEGMENTS
    y = ca.SX.sym('y', 3 * (2*n + 1))
    blocks = [y[3*i:3*i + 3] for i in range(2*n + 1)]
    pos = [blocks[pos_block(k)] for k in range(n + 1)]
    vel = [ca.SX.zeros(3)] + [blocks[vel_block(k)] for k in range(1, n + 1)]

    force = [calc_spring_force(pos[k], pos[k+1], vel[k], vel[k+1]) for k in range(n)]
    ydot = [ca.SX.zeros(3)] * (2*n + 1)
    for k in range(1, n + 1):
        above = force[k] if k < n else ca.SX.zeros(3)
        ydot[pos_block(k)] = vel[k]
        ydot[vel_block(k)] = G_EARTH - (force[k-1] - above) / MASS

    y0 = np.zeros(3 * (2*n + 1))
    for k in range(n + 1):
        y0[3*pos_block(k):3*pos_block(k) + 3] = [0.0, 0.0, -k * L_0]
        if k > 0:
            y0[3*vel_block(k):3*vel_block(k) + 3] = [0.0, 0.0, k * V0 / n]
    return y, ca.vertcat(*ydot), y0


def run_example():
    y, ydot, y0 = build_model()
    tfinal = 10.0           # must match Tether_04.jl's duration
    ncp    = 500            # number of communication points, must match Tether_04.jl's dt
    time = np.linspace(0.0, tfinal, ncp + 1)
    sim = ca.integrator('sim', 'cvodes', {'x': y, 'ode': ydot}, 0.0, time[1:],
                        {'abstol': 1.0e-6, 'reltol': 1.0e-6})
    y_sol = np.column_stack([y0, np.array(sim(x0=y0)['xf'])]).T

    # extract the z position and velocity of the lowest mass (mass SEGMENTS)
    pos_z_ix = 5 + (SEGMENTS - 1) * 6
    vel_z_ix = pos_z_ix + 3
    pos_z = y_sol[:, pos_z_ix]
    vel_z = y_sol[:, vel_z_ix]

    # saving the result for comparison with the Julia implementation
    os.makedirs("output", exist_ok=True)
    with open(os.path.join("output", "Tether_04_python.csv"), "w") as f:
        f.write("time,pos_z,vel_z\n")
        for t_i, pz_i, vz_i in zip(time, pos_z, vel_z):
            f.write(f"{t_i},{pz_i},{vz_i}\n")

    plt.gcf().canvas.manager.set_window_title("segmented tether")
    plt.ax1 = plt.subplot(111)
    plt.ax1.set_xlabel('time [s]')
    plt.plot(time, pos_z, color="green")
    plt.ax1.set_ylabel('pos_z [m]')
    plt.ax1.grid(True)
    plt.ax2 = plt.twinx()
    plt.ax2.set_ylabel('vel_z [m/s]')
    plt.plot(time, vel_z, color="red")

    if os.environ.get("TETHERS_BRIEF_PLOT") == "1":
        # show briefly and close automatically, e.g. when running the tests
        plt.pause(1)
        plt.close('all')
    else:
        plt.show()

if __name__ == '__main__':
    run_example()
