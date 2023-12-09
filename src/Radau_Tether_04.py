# -*- coding: utf-8 -*-
"""
Tutorial example simulating a 3D mass-spring system.
"""
import numpy as np
import pylab as plt
from assimulo.problem import Implicit_Problem #Imports the problem formulation from Assimulo
from assimulo.solvers import Radau5DAE #Imports the solver IDA from Assimulo

G_EARTH  = np.array([0.0, 0.0, -9.81]) # gravitational acceleration
C_SPRING = 50.0                        # spring constant [N/m] 
DAMPING  =  0.5                        # damping [Ns/m]
L_0      = 10.0                        # initial segment length [m]
SEGMENTS = 3
MASS     = 1.0                         # mass per tether of initial segment length [kg]
ZEROS  = np.array([0.0, 0.0, 0.0])
RESULT = np.zeros(SEGMENTS * 6 + 3).reshape((-1, 3))

# Example two:
# Three falling masses, connected with two springs with damping
# State vector y   = mass0.pos, mass1.pos, mass1.vel, mass2.pos, mass2.vel
# Derivative   yd  = mass0.vel, mass1.vel, mass1.acc, mass2.vel, mass2.acc
# Residual     res = (yd.mass0.vel), (y.mass1.vel - yd.mass1.vel), (yd.mass1.acc - G_EARTH),
#                                    (y.mass2.vel - yd.mass2.vel), (yd.mass2.acc - G_EARTH)    
def res1(t, y, yd):
    y1  = y.reshape((-1, 3)) # reshape the state vector such that we can access it per 3D-vector
    yd1 = yd.reshape((-1, 3))
    RESULT[0] = y1[0] # the velocity of mass0 shall be zero
    
    last_force = ZEROS # later this shall be the kite force
    
    for i in range(SEGMENTS-2, -1, -1):    # count down from segments-2 to zero
        # 1. calculate the force of the lowest spring (the spring next to the kite)   
        res_3 =  y1[2*i+4]  - yd1[2*i+3]  # the derivative of the position of mass1 must be equal to its velocity
        rel_vel1 = yd1[2*i+3] - yd1[2*i+1] # calcultate the relative velocity of mass2 with respect to mass 1 
        segment1 = y1[2*i+3] - y1[2*i+1]  # calculate the vector from mass1 to mass0
        force = C_SPRING * (np.linalg.norm(segment1) - L_0) * segment1 / np.linalg.norm(segment1) \
                         + DAMPING * rel_vel1 
                             
        # 2. apply it to the lowest mass (the mass next to the kite)   
        spring_forces = force - last_force    
        last_force = force                         
        acc1 = spring_forces / MASS  # create the vector of the spring acceleration    
        res_4 = yd1[2*i+4] - (G_EARTH - acc1) # the derivative of the velocity must be equal to the total acceleration  
        RESULT[2*i+3] = res_3 
        RESULT[2*i+4] = res_4                       

    # 3. calculate the force of the spring above    
    res_1 = y1[2]  - yd1[1]  # the derivative of the position of mass1 must be equal to its velocity
    rel_vel = yd1[1] - yd1[0] # calcultate the relative velocity of mass1 with respect to mass 0 
    segment = y1[1] - y1[0]  # calculate the vector from mass1 to mass0
    force = C_SPRING * (np.linalg.norm(segment) - L_0) * segment / np.linalg.norm(segment) \
                     + DAMPING * rel_vel    

    # 2. apply it to the next mass nearer to the winch
    spring_forces = force - last_force    
    acc = spring_forces / MASS  # create the vector of the spring acceleration
    res_2 = yd1[2] - (G_EARTH - acc) # the derivative of the velocity must be equal to the total acceleration
    
    RESULT[1] = res_1 
    RESULT[2] = res_2 
    return RESULT.flatten()

def run_example():
    # Set the initial conditons
    t0  = 0.0                   # Initial time
    pos, vel, acc = [], [], []
    x, y, z0 = 0.0, 0.0, 0.0
    dz = L_0
    for i in range (SEGMENTS + 1):
        z = z0 - i * dz
        if i == 0:
            pos.append(np.array([0.0, 0.0, z]))
        else:
            pos.append(np.array([x, y, z]))
        vel.append(np.array([0.0, 0.0, 0.0]))
        acc.append(np.array([0.0, 0.0, -9.81]))
    y0, yd0 = pos[0], vel[0]
    sw0 = []
    for i in range (SEGMENTS):    
        y0  = np.append(y0,  np.append(pos[i+1], vel[i+1])) # Initial state vector
        yd0 = np.append(yd0, np.append(vel[i+1], acc[i+1])) # Initial state vector derivative
        # array of booleans; true means the tether segment is loose (l < l_reelout)
        pos_ix = 6 * i  + 3 # position index of mass i + 1
        if i==0:
            last_pos_ix = 0
        else:
            last_pos_ix = pos_ix - 6        
        # calculate the norm of the vector from mass1 to mass0 minus the initial segment length
        event = (np.linalg.norm(y0[pos_ix:pos_ix + 3] - y0[last_pos_ix:last_pos_ix + 3] ) - L_0) <= 0
        sw0.append(event)       

    model = Implicit_Problem(res1, y0, yd0, t0) # Create an Assimulo problem
    model.name = 'Mass-Spring' # Specifies the name of problem (optional)

    sim = Radau5DAE(model)        # Create the IDA solver
        
    tfinal = 10.0           # Specify the final time
    ncp    = 4000            # Number of communcation points (number of return points)
    
    # Use the .simulate method to simulate and provide the final time and ncp (optional)    
    time, y, yd = sim.simulate(tfinal, ncp) 
    
    # plot the result
    pos_z1 = y[:,5]
    pos_z2 = y[:,5+6]   
    pos_z3 = y[:,5+2*6]  



    plt.ax1 = plt.subplot(111) 
    plt.ax1.set_xlabel('time [s]')
    plt.plot(time, pos_z1, color="green")
    plt.plot(time, pos_z2, color="blue")    
    plt.plot(time, pos_z3, color="black")    
    plt.ax1.set_ylabel('pos_z [m]')   
    if True:
        # vel_z = y[:,8]
        force_z = y[:,11]
        # rel_vel = yd[:,3:6] - yd[:,0:3]
        # vel_norm = np.sum(np.abs(rel_vel)**2, axis=-1)**(1./2)        
        plt.ax2 = plt.twinx()  
        plt.ax2.set_ylabel('vel_z [m/s]')   
        # plt.plot(time, vel_z, color="red")  
        plt.ax2.set_ylabel('force_z [N]')   
        plt.plot(time, force_z, color="red")      
        # plt.plot(time, vel_norm, color="blue")        
    plt.show()

if __name__ == '__main__':
    run_example()
