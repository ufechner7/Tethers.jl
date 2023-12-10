# -*- coding: utf-8 -*-
"""
Tutorial example simulating a 3D mass-spring system with a nonlinear spring (no spring forces
for l < l_0).
"""
import numpy as np
import pylab as plt
L_0 = 10.0 # initial segment length [m]

from tools import ExtProblem
from assimulo.solvers.sundials import IDA #Imports the solver IDA from Assimulo

G_EARTH  = np.array([0.0, 0.0, -9.81]) # gravitational acceleration
C_SPRING =  50                         # spring constant [N/m]; Dynema, 4mm: 200e3
DAMPING  =  0.5                        # damping [Ns/m]
MASS     = 1.0                         # mass per point-mass [kg]
V0       = 4.0                         # initial velocity

#Extend Assimulos problem definition
class ExtendedProblem(ExtProblem):
    # Set the initial conditions
    t0  = 0.0                   # Initial time
    pos_0 = np.array([0.0, 0.0, 0.0])     # Initial position of mass zero
    vel_0 = np.array([0.0, 0.0, 0.0])     # Initial velocity of mass zero
    pos_1 = np.array([0.0, 0.0,  -L_0])   # Initial position of mass one
    vel_1 = np.array([0.0, 0.0,  V0])    # Initial velocity of mass one
    acc_1 = np.array([0.0, 0.0, -9.81])   # Initial acceleration mass one
    y0  = np.append(pos_0, np.append(pos_1, vel_1)) # Initial state vector
    yd0 = np.append(vel_0, np.append(vel_1, acc_1)) # Initial state vector derivative
    sw0 = [vel_1[2] > 0] # array of booleans; true means the tether segment is loose (l < l_0)

    # Falling mass, attached to a spring
    # State vector y   = mass0.pos, mass1.pos, mass1.vel
    # Derivative   yd  = mass0.vel, mass1.vel, mass1.acc
    # Residual     res = (yd.mass0.vel), (y.mass1.vel - yd.mass1.vel), (yd.mass1.acc - G_EARTH)
    def res(self, t, y, yd, sw):
        res_0 = y[0:3]              # the velocity of mass0 shall be zero
        res_1 = y[6:9]  - yd[3:6]   # the derivative of the position of mass1 must be equal to its velocity
        rel_vel = yd[3:6] - yd[0:3] # calculate the relative velocity of mass1 with respect to mass 0
        segment = y[3:6] - y[0:3]   # calculate the vector from mass1 to mass0
        if not sw[0]:               # if the segment is not loose, calculate spring and damping force
            force = C_SPRING * (np.linalg.norm(segment) - L_0) * segment / np.linalg.norm(segment) \
                + DAMPING * rel_vel
        else:
            force = np.array([0.0, 0.0, 0.0])
        acc = force / MASS                # create the vector of the spring acceleration
        res_2 = yd[6:9] - (G_EARTH - acc) # the derivative of the velocity must be equal to the total acceleration
        return np.append(res_0, np.append(res_1, res_2))

    def state_events(self, t, y, yd, sw):
        """
        This is our function that keeps track of our events. When the sign
        of any of the events has changed, we have an event.
        """
        # calculate the norm of the vector from mass1 to mass0 minus the initial segment length
        event_0 = np.linalg.norm(y[3:6] - y[0:3]) - L_0
        return np.array([event_0])

def run_example():
    #Create an instance of the problem
    model = ExtendedProblem() #Create the problem
    model.name = 'Mass-Spring' # Specifies the name of problem (optional)
    sim = IDA(model) #Create the solver
    sim.verbosity = 30
    time, y, yd = sim.simulate(10.0, 500) #Simulate 10 seconds with 500 communications points

    # plot the result
    pos_z = y[:,5]
    vel_z = y[:,8]
    plt.ax1 = plt.subplot(111)
    plt.ax1.set_xlabel('time [s]')
    plt.plot(time, pos_z, color="green")
    plt.ax1.set_ylabel('pos_z [m]')
    plt.grid(True)
    plt.ax2 = plt.twinx()
    plt.ax2.set_ylabel('vel_z [m/s]')
    plt.plot(time, vel_z, color="red")
    plt.show()

if __name__ == '__main__':
    run_example()
