# -*- coding: utf-8 -*-
"""
Tutorial example simulating a falling mass, attached to a linear spring.
"""
import os
import numpy as np
import matplotlib.pyplot as plt

from assimulo.problem import Implicit_Problem # Imports the problem formulation from Assimulo
from assimulo.solvers.sundials import IDA # Imports the solver IDA from Assimulo

G_EARTH  = np.array([0.0, 0.0, -9.81]) # gravitational acceleration
C_SPRING =  50.0                       # spring constant [N/m]
DAMPING  =  0.5                        # damping [Ns/m]
MASS     = 1.0                         # mass of the point-mass [kg]
L_0      = 10.0                        # initial spring length [m]

# Extend Assimulos problem definition
class ExtendedProblem(Implicit_Problem):
    # Set the initial conditions
    t0  = 0.0                        # Initial time
    pos_0 = np.array([0.0, 0.0, -L_0]) # Initial position of the mass
    vel_0 = np.array([0.0, 0.0,  0.0]) # Initial velocity of the mass
    acc_0 = np.array([0.0, 0.0, -9.81]) # Initial acceleration of the mass
    y0  = np.append(pos_0, vel_0)    # Initial state vector
    yd0 = np.append(vel_0, acc_0)    # Initial state vector derivative

    # Falling mass, attached to a linear spring anchored at the origin
    # State vector y   = mass.pos, mass.vel
    # Derivative   yd  = mass.vel, mass.acc
    # Residual     res = (y.vel - yd.pos), (yd.vel - acc)
    def res(self, t, y, yd):
        res_0 = y[3:6] - yd[0:3]   # the derivative of the position must be equal to the velocity
        pos = y[0:3]
        vel = y[3:6]
        norm1 = np.linalg.norm(pos)
        unit_vector = -pos / norm1                     # direction from point mass to origin
        spring_vel = -np.dot(unit_vector, vel)
        spring_force = (C_SPRING * (norm1 - abs(L_0)) + DAMPING * spring_vel) * unit_vector
        acc = G_EARTH + spring_force / MASS
        res_1 = yd[3:6] - acc      # the derivative of the velocity must be equal to the total acceleration
        return np.append(res_0, res_1)

def run_example():
    # Create an instance of the problem
    model = ExtendedProblem()   # Create the problem
    model.name = 'Falling mass, linear spring' # Specifies the name of problem (optional)
    sim = IDA(model)            # Create the solver
    sim.verbosity = 30
    # let the solver pick its own steps, then resample onto the same time grid
    # as the Julia implementation (0:0.02:10) -- requesting communication points
    # directly can make Assimulo silently drop a point close to t_final
    raw_time, raw_y, raw_yd = sim.simulate(10.0)

    time = np.linspace(0.0, 10.0, 501)
    pos_z = np.interp(time, raw_time, raw_y[:,2])
    vel_z = np.interp(time, raw_time, raw_y[:,5])

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
