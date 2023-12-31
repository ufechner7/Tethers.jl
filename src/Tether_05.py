# -*- coding: utf-8 -*-
"""
Tutorial example simulating a 3D mass-spring system with a nonlinear spring (no spring forces
for l < l_0). It uses five tether segments. The force coupling is now implemented
correctly.
"""
import numpy as np
import pylab as plt
import math
import time

from assimulo.solvers.sundials import IDA     # Imports the solver IDA from Assimulo
from assimulo.problem import Implicit_Problem # Imports the problem formulation from Assimulo

G_EARTH  = np.array([0.0, 0.0, -9.81]) # gravitational acceleration
C_SPRING = 50.0                        # spring constant
DAMPING  =  0.5                        # damping [Ns/m]
L0      =  5.0                         # initial segment length     [m]
ALPHA0   = math.pi/10                  # initial tether angle     [rad]
SEGMENTS = 5
MASS     = 0.5                         # mass per tether particle  [kg]   
DURATION = 30                          # duration of the simulation [s]
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
        l0 = -i*L0
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
    print(y0)
    print(yd0)

    def res(self, t, y, yd):  
        y1  = y.reshape((-1, 3)) # reshape the state vector such that we can access it per 3D-vector
        yd1 = yd.reshape((-1, 3))
        length = L0
        c_spring = C_SPRING        
        damping  = DAMPING
        RESULT[0] = y1[0] # the velocity of mass0 shall be zero        
        last_force = ZEROS # later this shall be the kite force
        
        for i in range(SEGMENTS-2, -1, -1):    # count down from segments-2 to zero
            # 1. calculate the force of the lowest spring (the spring next to the kite)   
            res_3   =  y1[2*i+4] - yd1[2*i+3]  # the derivative of the position of mass1 must be equal to its velocity
            rel_vel = yd1[2*i+3] - yd1[2*i+1]  # calculate the relative velocity of mass2 with respect to mass 1 
            segment = y1[2*i+3]  - y1[2*i+1]   # calculate the vector from mass1 to mass0
            if np.linalg.norm(segment) > L0:               # if the segment is not loose, calculate spring and damping force
                c_spring = C_SPRING
            else:
                c_spring = 0.0
            force = c_spring * (np.linalg.norm(segment) - L0) * segment / np.linalg.norm(segment) \
                    + damping * rel_vel                                                
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
        norm = math.sqrt(segment[0]**2 + segment[1]**2 + segment[2]**2)  
        force = c_spring * (norm - length) * segment / norm + damping * rel_vel    
    
        # 2. apply it to the next mass nearer to the winch
        spring_forces = force - last_force    
        mass = MASS          
        acc = spring_forces / mass  # create the vector of the spring acceleration
        res_2 = yd1[2] - (G_EARTH - acc) # the derivative of the velocity must be equal to the total acceleration
        
        RESULT[1] = res_1 
        RESULT[2] = res_2 
        return RESULT.flatten()
    
def plot2d(fig, y, reltime, segments, line, sc, txt):
    index = round(reltime*50)
    # print("index: ", index)
    x, z = np.zeros(segments+1), np.zeros(segments+1)
    for i in range(segments):
        x[i+1] = y[index, 3+6*i]
        z[i+1] = y[index, 5+6*i]
    z_max = np.max(z)
    if line is None:
        line, = plt.plot(x,z, linewidth="1")
        sc  = plt.scatter(x, z, s=15, color="red")
        plt.pause(0.01)
        txt = plt.annotate("t="+str(reltime)+" s",  
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

def play(duration, y):
    dt = 0.15
    plt.ioff()
    fig = plt.figure()
    plt.ylim(-SEGMENTS*L0-10, 0.5)
    plt.xlim(-SEGMENTS*L0/2, SEGMENTS*L0/2)
    plt.grid(True, color="grey", linestyle="dotted")
    line, sc, txt = None, None, None
    for t in np.linspace(0, duration, num=round(duration/dt)+1):
        line, sc, txt = plot2d(fig, y, t, SEGMENTS, line, sc, txt)
        time.sleep(dt/2)
    plt.show()
   
def run_example():  
    # Create an instance of the problem 
    model = ExtendedProblem()  # Create the problem 
    model.name = 'Mass-Spring' # Specifies the name of problem (optional)   
    
    sim = IDA(model) # Create the solver 
    sim.verbosity = 30 
    sim.atol = 1.0e-6
    sim.rtol = 1.0e-6
    
    time, y, yd = sim.simulate(DURATION, round(DURATION*50)+2) # 50 communications points per second 
    play(DURATION, y)   
    
if __name__ == '__main__':
    run_example()
