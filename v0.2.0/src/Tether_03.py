# -*- coding: utf-8 -*-
"""
Tutorial example simulating a 3D mass-spring system with a nonlinear spring (no spring forces
for l < l_0).
"""
import numpy as np
import pylab as plt

from assimulo.problem import Implicit_Problem # Imports the problem formulation from Assimulo
from assimulo.solvers.sundials import IDA # Imports the solver IDA from Assimulo

G_EARTH  = np.array([0.0, 0.0, -9.81]) # gravitational acceleration
C_SPRING =  50                         # spring constant [N/m]; Dynema, 4mm: 200e3
DAMPING  =  0.5                        # damping [Ns/m]
MASS     = 1.0                         # mass per point-mass [kg]
L_0      = 10.0                        # initial segment length [m]
V0       = 4.0                         # initial velocity

# Extend Assimulos problem definition
class ExtendedProblem(Implicit_Problem):
    # Set the initial conditions
    t0  = 0.0                   # Initial time
    pos_0 = np.array([0.0, 0.0, 0.0])     # Initial position of mass zero
    vel_0 = np.array([0.0, 0.0, 0.0])     # Initial velocity of mass zero
    pos_1 = np.array([0.0, 0.0,  -L_0])   # Initial position of mass one
    vel_1 = np.array([0.0, 0.0,  V0])     # Initial velocity of mass one
    acc_1 = np.array([0.0, 0.0, -9.81])   # Initial acceleration mass one
    y0  = np.append(pos_0, np.append(pos_1, vel_1)) # Initial state vector
    yd0 = np.append(vel_0, np.append(vel_1, acc_1)) # Initial state vector derivative

    # Falling mass, attached to a spring
    # State vector y   = mass0.pos, mass1.pos, mass1.vel
    # Derivative   yd  = mass0.vel, mass1.vel, mass1.acc
    # Residual     res = (yd.mass0.vel), (y.mass1.vel - yd.mass1.vel), (yd.mass1.acc - G_EARTH)
    def res(self, t, y, yd):
        res_0 = y[0:3]              # the velocity of mass0 shall be zero
        res_1 = y[6:9]  - yd[3:6]   # the derivative of the position of mass1 must be equal to its velocity
        rel_vel = yd[3:6] - yd[0:3] # calculate the relative velocity of mass1 with respect to mass 0
        segment = y[3:6] - y[0:3]   # calculate the vector from mass1 to mass0
        if np.linalg.norm(segment) > L_0:               # if the segment is not loose, calculate spring and damping force
            c_spring = C_SPRING
        else:
            c_spring = 0.0
        force = c_spring * (np.linalg.norm(segment) - L_0) * segment / np.linalg.norm(segment) \
                + DAMPING * rel_vel
        acc = force / MASS                # create the vector of the spring acceleration
        res_2 = yd[6:9] - (G_EARTH - acc) # the derivative of the velocity must be equal to the total acceleration
        return np.append(res_0, np.append(res_1, res_2))

def run_example():
    # Create an instance of the problem
    model = ExtendedProblem()  # Create the problem
    model.name = 'Mass-Spring' # Specifies the name of problem (optional)
    sim = IDA(model)           # Create the solver
    sim.verbosity = 30
    time, y, yd = sim.simulate(10.0, 500) #Simulate 10 seconds with 500 communications points
    print(len(time))

    # plot the result
    pos_z = y[:,5]
    vel_z = y[:,8]
    C_SPRINGS = np.zeros(len(time))
    for i in range(len(time)):
        yi=y[i]
        segment = yi[3:6] - yi[0:3]
        if np.linalg.norm(segment) > L_0:               # if the segment is not loose, calculate spring and damping force
            c_spring = C_SPRING
        else:
            c_spring = 0.0
        C_SPRINGS[i] = c_spring
    plt.ax1 = plt.subplot(111)
    plt.ax1.set_xlabel('time [s]')
    plt.plot(time, pos_z, color="green")
    plt.plot(time, -np.ones(len(time)) * L_0 + 0.005 * C_SPRINGS, color="grey", label="c_spring")
    plt.ax1.set_ylabel('pos_z [m]')
    plt.grid(True)
    plt.ax2 = plt.twinx()
    plt.ax2.set_ylabel('vel_z [m/s]')
    plt.plot(time, vel_z, color="red")
    plt.show()

if __name__ == '__main__':
    run_example()
