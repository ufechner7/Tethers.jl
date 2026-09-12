# -*- coding: utf-8 -*-
"""
Tutorial example simulating a falling mass, attached to a linear spring.

The model is written as a CasADi expression graph and integrated with SUNDIALS' CVODES,
which CasADi supplies with the exact Jacobian of the model.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import casadi as ca

G_EARTH  = np.array([0.0, 0.0, -9.81]) # gravitational acceleration
C_SPRING =  50.0                       # spring constant [N/m]
DAMPING  =  0.5                        # damping [Ns/m]
MASS     = 1.0                         # mass of the point-mass [kg]
L_0      = 10.0                        # initial spring length [m]

# Falling mass, attached to a linear spring anchored at the origin
# State vector y = mass.pos, mass.vel
def build_model():
    """ The mass and its spring as one CasADi expression graph. Returns the state vector
        `y`, its derivative `ydot` and the initial state `y0`. """
    y = ca.SX.sym('y', 6)
    pos, vel = y[0:3], y[3:6]
    norm1 = ca.norm_2(pos)
    unit_vector = -pos / norm1                     # direction from point mass to origin
    spring_vel = -ca.dot(unit_vector, vel)
    spring_force = (C_SPRING * (norm1 - abs(L_0)) + DAMPING * spring_vel) * unit_vector
    acc = G_EARTH + spring_force / MASS
    y0 = np.array([0.0, 0.0, -L_0, 0.0, 0.0, 0.0])   # pos, vel
    return y, ca.vertcat(vel, acc), y0

def run_example():
    y, ydot, y0 = build_model()
    time = np.linspace(0.0, 10.0, 501)
    sim = ca.integrator('sim', 'cvodes', {'x': y, 'ode': ydot}, 0.0, time[1:],
                        {'abstol': 1.0e-6, 'reltol': 1.0e-6})
    y_sol = np.column_stack([y0, np.array(sim(x0=y0)['xf'])]).T
    pos_z = y_sol[:, 2]
    vel_z = y_sol[:, 5]

    # saving the result for comparison with the Julia implementation
    os.makedirs("output", exist_ok=True)
    with open(os.path.join("output", "Tether_02_python.csv"), "w") as f:
        f.write("time,pos_z,vel_z\n")
        for t_i, pz_i, vz_i in zip(time, pos_z, vel_z):
            f.write(f"{t_i},{pz_i},{vz_i}\n")

    plt.gcf().canvas.manager.set_window_title("falling mass, linear spring")
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
