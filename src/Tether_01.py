# -*- coding: utf-8 -*-
"""
Tutorial example showing how to use the implicit solver RADAU. 
It simulates a falling mass.
"""
import numpy as np
import pylab as plt
from assimulo.problem import Implicit_Problem # Imports the problem formulation from Assimulo
from assimulo.solvers import Radau5DAE        # Imports the solver RADAU from Assimulo

G_EARTH  = np.array([0.0, 0.0, -9.81]) # gravitational acceleration
  
# Example one:
# Falling mass.
# State vector y   = mass.pos, mass.vel
# Derivative   yd  = mass.vel, mass.acc
# Residual     res = (y.vel - yd.vel), (yd.acc - G_EARTH)     
def res1(t, y, yd):
    res_0 = y[3:6]  - yd[0:3]
    res_1 = yd[3:6] - G_EARTH 
    return np.append(res_0, res_1)

def run_example():
    # Set the initial conditons
    t0  = 0.0                   # Initial time
    vel_0 = [0.0, 0.0, 50.0]    # Initial velocity
    pos_0 = [0.0, 0.0,  0.0]    # Initial position
    acc_0 = [0.0, 0.0, -9.81]   # Initial acceleration
    y0  = pos_0 + vel_0         # Initial pos, vel
    yd0 = vel_0 + acc_0         # Initial vel, acc

    model = Implicit_Problem(res1, y0, yd0, t0) # Create an Assimulo problem
    model.name = 'Falling mass' # Specifies the name of problem (optional)

    sim = Radau5DAE(model)      # Create the IDA solver
        
    tfinal = 10.0           # Specify the final time
    ncp    = 500            # Number of communcation points (number of return points)

    # Use the .simulate method to simulate and provide the final time and ncp (optional)    
    time, y, yd = sim.simulate(tfinal, ncp) 
    
    # plot the result
    pos_z = y[:,2]
    vel_z = y[:,5]
    plt.ax1 = plt.subplot(111) 
    plt.ax1.set_xlabel('time [s]')
    plt.plot(time, pos_z, color="green")
    plt.ax1.set_ylabel('pos_z [m]')  
    plt.ax1.grid(True) 
    plt.ax2 = plt.twinx()  
    plt.ax2.set_ylabel('vel_z [m/s]')   
    plt.plot(time, vel_z, color="red")    
    plt.show()

if __name__ == '__main__':
    run_example()
    