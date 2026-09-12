# -*- coding: utf-8 -*-
"""
Tutorial example showing how to use an implicit solver. It simulates a falling mass.

The model is written as a CasADi expression graph and integrated with SUNDIALS' CVODES,
which CasADi supplies with the exact Jacobian of the model.
"""
import os
import numpy as np
import pylab as plt
import casadi as ca

G_EARTH  = np.array([0.0, 0.0, -9.81]) # gravitational acceleration

# Example one: Falling mass
# State vector y = mass.pos, mass.vel
def build_model():
    """ The falling mass as one CasADi expression graph. Returns the state vector `y`,
        its derivative `ydot` and the initial state `y0`. """
    y = ca.SX.sym('y', 6)
    ydot = ca.vertcat(y[3:6], ca.DM(G_EARTH))
    y0 = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 50.0])   # pos, vel
    return y, ydot, y0

def run_example():
    y, ydot, y0 = build_model()
    tfinal = 10.0           # Specify the final time
    ncp    = 500            # Number of communication points (number of return points)
    time = np.linspace(0.0, tfinal, ncp + 1)
    sim = ca.integrator('sim', 'cvodes', {'x': y, 'ode': ydot}, 0.0, time[1:],
                        {'abstol': 1.0e-6, 'reltol': 1.0e-6})
    y_sol = np.column_stack([y0, np.array(sim(x0=y0)['xf'])]).T

    # plot the result
    pos_z = y_sol[:, 2]
    vel_z = y_sol[:, 5]

    # saving the result for comparison with the Julia implementation
    os.makedirs("output", exist_ok=True)
    with open(os.path.join("output", "Tether_01_python.csv"), "w") as f:
        f.write("time,pos_z,vel_z\n")
        for t_i, pz_i, vz_i in zip(time, pos_z, vel_z):
            f.write(f"{t_i},{pz_i},{vz_i}\n")

    plt.gcf().canvas.manager.set_window_title("falling mass")
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
    