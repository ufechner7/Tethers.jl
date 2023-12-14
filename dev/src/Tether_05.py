# -*- coding: utf-8 -*-
"""
Tutorial example simulating a 3D mass-spring system with a nonlinear spring (no spring forces
for l < l_0). It uses five tether segments. The force coupling is now implemented
correctly.
"""
import numpy as np
import pylab as plt
import math

from assimulo.solvers.sundials import IDA     # Imports the solver IDA from Assimulo
from assimulo.problem import Implicit_Problem # Imports the problem formulation from Assimulo

G_EARTH  = np.array([0.0, 0.0, -9.81]) # gravitational acceleration
C_SPRING = 50.0                        # spring constant
DAMPING  =  0.5                        # damping [Ns/m]
L_0      =  5.0                        # initial segment length    [m]
ALPHA0   = math.pi/8                   # initial tether angle    [rad]
SEGMENTS = 10
MASS     = 0.5                         # mass per tether particle [kg]      
ZEROS  = np.array([0.0, 0.0, 0.0])
RESULT = np.zeros(SEGMENTS * 6 + 3).reshape((-1, 3))
NONLINEAR = True

# State vector y   = mass0.pos, mass1.pos, mass1.vel
# Derivative   yd  = mass0.vel, mass1.vel, mass1.acc
# Residual     res = (yd.mass0.vel), (y.mass1.vel - yd.mass1.vel), (yd.mass1.acc - G_EARTH)     

# Extend Assimulos problem definition
class ExtendedProblem(Implicit_Problem):
    # Set the initial conditions
    t0  = 0.0                   # Initial time
    pos, vel, acc = [], [], []
    x, y, z0 = 0.0, 0.0, 0.0
    for i in range (SEGMENTS + 1):
        l0 = -i*L_0
        pos.append(np.array([math.sin(ALPHA0) * l0, 0.0, math.cos(ALPHA0) * l0]))            
        vel.append(np.array([0.0, 0.0, 0.0]))
        if i == 0:
            acc.append(np.array([0.0, 0.0, 0.0]))
        else:
            acc.append(np.array([0.0, 0.0, -9.81]))
    y0, yd0 = pos[0], vel[0]
    for i in range (SEGMENTS):    
        y0  = np.append(y0,  np.append(pos[i+1], vel[i+1])) # Initial state vector
        yd0 = np.append(yd0, np.append(vel[i+1], acc[i+1])) # Initial state vector derivative
        pos_ix = 6 * i  + 3 # position index of mass i + 1
        if i==0:
            last_pos_ix = 0
        else:
            last_pos_ix = pos_ix - 6            
    print(y0)
    print(yd0)

    def res(self, t, y, yd, sw):  
        y1  = y.reshape((-1, 3)) # reshape the state vector such that we can access it per 3D-vector
        yd1 = yd.reshape((-1, 3))
        length = L_0
        c_spring = C_SPRING * L_0 / length        
        damping = DAMPING * length / L_0
        RESULT[0] = y1[0] # the velocity of mass0 shall be zero        
        last_force = ZEROS # later this shall be the kite force
        
        for i in range(SEGMENTS-2, -1, -1):    # count down from segments-2 to zero
            # 1. calculate the force of the lowest spring (the spring next to the kite)   
            res_3   =  y1[2*i+4] - yd1[2*i+3]  # the derivative of the position of mass1 must be equal to its velocity
            rel_vel = yd1[2*i+3] - yd1[2*i+1]  # calculate the relative velocity of mass2 with respect to mass 1 
            segment = y1[2*i+3]  - y1[2*i+1]   # calculate the vector from mass1 to mass0
            if not sw[SEGMENTS - 1 - i] or not NONLINEAR:               # if the segment is not loose, calculate spring and damping force
                norm = math.sqrt(segment[0]**2 + segment[1]**2 + segment[2]**2)                
                force = c_spring * (norm - length) * segment / norm + damping * rel_vel 
            else:
                force = damping * rel_vel                                                    
            # 2. apply it to the lowest mass (the mass next to the kite)   
            spring_forces = force - last_force    
            last_force = force       
            mass = MASS               
            acc1 = spring_forces / mass  # create the vector of the spring acceleration    
            res_4 = yd1[2*i+4] - (G_EARTH - acc1) # the derivative of the velocity must be equal to the total acceleration  
            RESULT[2*i+3] = res_3 
            RESULT[2*i+4] = res_4                       
    
        # 3. calculate the force of the spring above    
        res_1   = y1[2]  - yd1[1] # the derivative of the position of mass1 must be equal to its velocity
        rel_vel = yd1[1] - yd1[0] # calculate the relative velocity of mass1 with respect to mass 0 
        segment = y1[1]  - y1[0]  # calculate the vector from mass1 to mass0
        if not sw[0] or not NONLINEAR:         
            norm = math.sqrt(segment[0]**2 + segment[1]**2 + segment[2]**2)  
            force = c_spring * (norm - length) * segment / norm + damping * rel_vel    
        else:
            force = ZEROS
    
        # 2. apply it to the next mass nearer to the winch
        spring_forces = force - last_force    
        mass = MASS          
        acc = spring_forces / mass  # create the vector of the spring acceleration
        res_2 = yd1[2] - (G_EARTH - acc) # the derivative of the velocity must be equal to the total acceleration
        
        RESULT[1] = res_1 
        RESULT[2] = res_2 
        return RESULT.flatten()
   
def run_example():  
    # Create an instance of the problem 
    model = ExtendedProblem()  # Create the problem 
    model.name = 'Mass-Spring' # Specifies the name of problem (optional)   
    
    sim = IDA(model) # Create the solver 
    sim.verbosity = 30 
    sim.atol = 1.0e-4
    sim.rtol = 1.0e-5
    #sim.continuous_output = True  
    
    time, y, yd = sim.simulate(20.0, 400) #Simulate 10 seconds with 1000 communications points     
    
    # plot the result
    pos_z1 = y[:,5]
    vel_z = y[:,8]
#    force_z = y[:,11]
    # rel_vel = yd[:,3:6] - yd[:,0:3]
    # vel_norm = np.sum(np.abs(rel_vel)**2, axis=-1)**(1./2)
    plt.ax1 = plt.subplot(111) 
    plt.ax1.set_xlabel('time [s]')
    plt.plot(time, pos_z1, color="green")
    if SEGMENTS > 1:    
        pos_z2 = y[:,5+6]       
        plt.plot(time, pos_z2, color="blue")    
    if SEGMENTS > 2:    
        pos_z3 = y[:,5+12]       
        plt.plot(time, pos_z3, color="yellow")      
    if SEGMENTS > 3:    
        pos_z4 = y[:,5+18]       
        plt.plot(time, pos_z4, color="grey")  
    if SEGMENTS > 4:    
        pos_z5 = y[:,5+24]     
        vel_z = y[:,8+24]
        plt.plot(time, pos_z5, color="black")            
    plt.ax1.set_ylabel('pos_z [m]')   
    plt.ax2 = plt.twinx()  
    plt.ax2.set_ylabel('vel_z [m/s]')   
    plt.plot(time, vel_z, color="red")  
#    plt.ax2.set_ylabel('force_z [N]')   
#    plt.plot(time, force_z, color="red")      
    # plt.plot(time, vel_norm, color="blue")  
    plt.grid(True)      
    plt.show()

if __name__ == '__main__':
    model = ExtendedProblem()  # Create the problem 
    # run_example()
