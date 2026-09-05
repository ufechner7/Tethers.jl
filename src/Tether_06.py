# -*- coding: utf-8 -*-
"""
Tutorial example simulating a 3D mass-spring system with a nonlinear spring (no spring forces
for l < l_0). It uses five tether segments. The force coupling is now implemented
correctly.
"""
import numpy as np
import matplotlib.pyplot as plt
import math
import time

from assimulo.solvers.sundials import IDA     # Imports the solver IDA from Assimulo
from assimulo.problem import Implicit_Problem # Imports the problem formulation from Assimulo

G_EARTH  = np.array([0.0, 0.0, -9.81]) # gravitational acceleration
C_SPRING =  614600.0                   # spring constant
DAMPING  =  473*0.0045                 # damping [Ns/m]
L0      =  50.0                        # initial segment length     [m]
ALPHA0   = math.pi/10                  # initial tether angle     [rad]
SEGMENTS = 5 
DURATION = 10                          # duration of the simulation [s]

V_RO = 2.0                             # reel-out speed                  [m/s]
D_TETHER = 4.0                         # tether diameter                  [mm]
RHO_TETHER = 724.0                     # density of Dyneema            [kg/m³] 

ZEROS  = np.array([0.0, 0.0, 0.0])
RESULT = np.zeros(SEGMENTS * 6 + 3).reshape((-1, 3))
mass_per_meter = RHO_TETHER * math.pi * (D_TETHER/2000.0)**2

def calc_spring_force(pos1, pos2, vel1, vel2, l_seg, c_spring0, damping):
    """ Spring and damping force of the segment between the masses at pos1 and pos2.
        The result points from pos1 towards pos2. Spring and damper act in parallel
        along the segment, therefore the damping force uses the component of the
        relative velocity along the segment and not the full 3D vector. """
    segment     = pos2 - pos1
    norm        = np.linalg.norm(segment)
    unit_vector = segment / norm
    c_spring    = c_spring0 if norm > l_seg else 0.0 # a loose segment cannot push
    spring_vel  = np.dot(vel2 - vel1, unit_vector)   # rate of change of the segment length
    return (c_spring * (norm - l_seg) + damping * spring_vel) * unit_vector

def calc_particle_mass(i, m_tether_particle):
    """ Mass of particle i. The last particle is attached to one segment only and
        therefore carries half of the mass of an inner particle. """
    if i == SEGMENTS:
        return 0.5 * m_tether_particle
    return m_tether_particle

# State vector y   = mass0.pos, mass1.pos, mass1.vel
# Derivative   yd  = mass0.vel, mass1.vel, mass1.acc
# Residual     res = (yd.mass0.vel), (y.mass1.vel - yd.mass1.vel), (yd.mass1.acc - G_EARTH)     

# Extend Assimulos problem definition
class ExtendedProblem(Implicit_Problem):
    # Set the initial conditions
    t0  = 0.0                   # Initial time
    pos, vel, acc = [], [], []
    for i in range (SEGMENTS + 1):
        l0 = -i*L0/SEGMENTS
        pos.append(np.array([math.sin(ALPHA0) * l0, 0.0, math.cos(ALPHA0) * l0]))            
        vel.append(np.array([0.0, 0.0, 0.0]))
        if i == 0:
            acc.append(np.array([0.0, 0.0, 0.0]))
        else:
            acc.append(G_EARTH)
    y0, yd0 = pos[0], vel[0]
    for i in range (SEGMENTS):    
        y0  = np.append(y0,  np.append(pos[i+1], vel[i+1])) # Initial state vector
        yd0 = np.append(yd0, np.append(vel[i+1], acc[i+1])) # Initial state vector derivative          

    def res(self, t, y, yd):  
        y1  = y.reshape((-1, 3)) # reshape the state vector such that we can access it per 3D-vector
        yd1 = yd.reshape((-1, 3))
        l_seg = (L0 + V_RO*t) / SEGMENTS   # unstretched length of one segment
        c_spring1 = C_SPRING / l_seg
        damping  = DAMPING / l_seg
        m_tether_particle = mass_per_meter * l_seg
        # mass0 is fixed: constrain its position. Its derivative yd1[0] is not
        # determined by any residual and stays at the initial guess of zero.
        RESULT[0] = y1[0]
        last_force = ZEROS # later this shall be the kite force
        
        for i in range(SEGMENTS-2, -1, -1):    # count down from segments-2 to zero
            # 1. calculate the force of the lowest spring (the spring next to the kite)   
            res_3   =  y1[2*i+4] - yd1[2*i+3]  # the derivative of the position of mass1 must be equal to its velocity
            # the force of the segment between the masses i+1 and i+2
            force = calc_spring_force(y1[2*i+1], y1[2*i+3], yd1[2*i+1], yd1[2*i+3],
                                      l_seg, c_spring1, damping)
            # 2. apply it to the lowest mass (the mass next to the kite)   
            spring_forces = force - last_force    
            last_force = force       
            mass = calc_particle_mass(i+2, m_tether_particle)
            acc1 = spring_forces / mass  # create the vector of the spring acceleration    
            res_4 = yd1[2*i+4] - (G_EARTH - acc1) # the derivative of the velocity must be equal to the total acceleration  
            RESULT[2*i+3] = res_3 
            RESULT[2*i+4] = res_4                       
    
        # 3. calculate the force of the spring above    
        res_1   = y1[2]  - yd1[1] # the derivative of the position of mass1 must be equal to its velocity
        force = calc_spring_force(y1[0], y1[1], yd1[0], yd1[1], l_seg, c_spring1, damping)

        # 2. apply it to the next mass nearer to the winch
        spring_forces = force - last_force    
        mass = calc_particle_mass(1, m_tether_particle)
        acc = spring_forces / mass  # create the vector of the spring acceleration
        res_2 = yd1[2] - (G_EARTH - acc) # the derivative of the velocity must be equal to the total acceleration
        
        RESULT[1] = res_1 
        RESULT[2] = res_2 
        return RESULT.flatten()
    
def plot2d(fig, t_sol, y, reltime, segments, line, sc, txt):
    index = min(np.searchsorted(t_sol, reltime), len(t_sol) - 1)
    x, z = np.zeros(segments+1), np.zeros(segments+1)
    for i in range(segments):
        x[i+1] = y[index, 3+6*i]
        z[i+1] = y[index, 5+6*i]
    if line is None:
        z_max = np.max(z)
        line, = plt.plot(x, z, linewidth=1)
        sc  = plt.scatter(x, z, s=15, color="red")
        plt.pause(0.01)
        txt = plt.annotate("t="+str(round(reltime,1))+" s",  
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

def play(duration, t_sol, y):
    dt = 0.151
    plt.ioff()
    fig = plt.figure()
    plt.ylim(-1.2*(L0+V_RO*duration), 0.5)
    plt.xlim(-L0/2, L0/2)
    plt.grid(True, color="grey", linestyle="dotted")
    plt.tight_layout(rect=(0, 0, 0.98, 0.98))
    line, sc, txt = None, None, None
    for t in np.linspace(0, duration, num=round(duration/dt)):
        line, sc, txt = plot2d(fig, t_sol, y, t, SEGMENTS, line, sc, txt)
        time.sleep(dt/2)
    plt.show(block=True)
   
def run_example():  
    # Create an instance of the problem 
    model = ExtendedProblem()  # Create the problem 
    model.name = 'Mass-Spring' # Specifies the name of problem (optional)   
    
    sim = IDA(model) # Create the solver 
    sim.verbosity = 0 
    sim.atol = 1.0e-6
    sim.rtol = 1.0e-6
    sim.linear_solver="SPGMR"

    t_sol, y, yd = sim.simulate(DURATION, round(DURATION*50)) # 50 communication points per second
    play(DURATION, t_sol, y)
    return
    
if __name__ == '__main__':
    run_example()
