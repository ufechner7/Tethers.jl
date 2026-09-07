# -*- coding: utf-8 -*-
"""
Tutorial example simulating a 3D mass-spring system with a nonlinear spring (no spring forces
for l < l_0). It uses five tether segments. The force coupling is now implemented
correctly.
"""
import numpy as np
import matplotlib.pyplot as plt
import math
import os
import time

from assimulo.solvers.sundials import IDA     # Imports the solver IDA from Assimulo
from assimulo.problem import Implicit_Problem # Imports the problem formulation from Assimulo

G_EARTH  = np.array([0.0, 0.0, -9.81]) # gravitational acceleration
C_SPRING = 50.0                        # spring constant
DAMPING  =  0.5                        # damping [Ns/m]
L0      =  5.0                         # initial segment length     [m]
V0       =  2.0                        # initial velocity of the lowest mass [m/s]
ALPHA0   = math.pi/10                  # initial tether angle     [rad]
SEGMENTS = 5
MASS     = 0.5                         # mass per tether particle  [kg]   
DURATION = 10                          # duration of the simulation [s]
ZEROS  = np.array([0.0, 0.0, 0.0])
RESULT = np.zeros(SEGMENTS * 6 + 3).reshape((-1, 3))
NONLINEAR = True                       # if True, a loose segment exerts no spring force

def calc_spring_constant(norm):
    """ Spring constant of one segment. A loose segment (norm <= L0) cannot push,
        therefore its spring constant is zero if NONLINEAR is True. """
    if NONLINEAR and norm <= L0:
        return 0.0
    return C_SPRING

def calc_spring_force(pos1, pos2, vel1, vel2):
    """ Spring and damping force of the segment between the masses at pos1 and pos2.
        The result points from pos1 towards pos2. Spring and damper act in parallel
        along the segment, therefore the damping force uses the component of the
        relative velocity along the segment and not the full 3D vector. """
    segment     = pos2 - pos1
    norm        = np.linalg.norm(segment)
    unit_vector = segment / norm
    spring_vel  = np.dot(vel2 - vel1, unit_vector) # rate of change of the segment length
    return (calc_spring_constant(norm) * (norm - L0) + DAMPING * spring_vel) * unit_vector

def calc_particle_mass(i):
    """ Mass of particle i. The last particle is attached to one segment only and
        therefore carries half of the mass of an inner particle. """
    if i == SEGMENTS:
        return 0.5 * MASS
    return MASS

def calc_acc(pos, vel):
    """ Acceleration of each particle for the given state. Particle zero is fixed.
        Used to calculate initial conditions that are physically consistent. """
    force = [calc_spring_force(pos[i], pos[i+1], vel[i], vel[i+1]) for i in range(SEGMENTS)]
    acc = [ZEROS]
    for i in range(1, SEGMENTS + 1):
        force_below = force[i] if i < SEGMENTS else ZEROS # no segment below the last particle
        acc.append(G_EARTH + (force_below - force[i-1]) / calc_particle_mass(i))
    return acc

# State vector y   = mass0.pos, mass1.pos, mass1.vel
# Derivative   yd  = mass0.vel, mass1.vel, mass1.acc
# Residual     res = (yd.mass0.vel), (y.mass1.vel - yd.mass1.vel), (yd.mass1.acc - G_EARTH)     

# Extend Assimulos problem definition
class ExtendedProblem(Implicit_Problem):
    # Set the initial conditions
    t0  = 0.0                   # Initial time
    pos, vel = [], []
    for i in range (SEGMENTS + 1):
        l0 = -i*L0
        v0 =  i*V0/SEGMENTS
        pos.append(np.array([math.sin(ALPHA0) * l0, 0.0, math.cos(ALPHA0) * l0]))            
        vel.append(np.array([math.sin(ALPHA0) * v0, 0.0, math.cos(ALPHA0) * v0]))
    acc = calc_acc(pos, vel) # consistent initial accelerations for the given pos and vel
    y0, yd0 = pos[0], vel[0]
    for i in range (SEGMENTS):    
        y0  = np.append(y0,  np.append(pos[i+1], vel[i+1])) # Initial state vector
        yd0 = np.append(yd0, np.append(vel[i+1], acc[i+1])) # Initial state vector derivative          
    print(y0)
    print(yd0)

    def res(self, t, y, yd):  
        y1  = y.reshape((-1, 3)) # reshape the state vector such that we can access it per 3D-vector
        yd1 = yd.reshape((-1, 3))
        # mass0 is fixed: constrain its position. Its derivative yd1[0] is not
        # determined by any residual and stays at the initial guess of zero.
        RESULT[0] = y1[0]
        last_force = ZEROS # later this shall be the kite force
        
        for i in range(SEGMENTS-2, -1, -1):    # count down from segments-2 to zero
            # 1. calculate the force of the lowest spring (the spring next to the kite)   
            res_3   =  y1[2*i+4] - yd1[2*i+3]  # the derivative of the position of mass1 must be equal to its velocity
            # the force of the segment between the masses i+1 and i+2
            force = calc_spring_force(y1[2*i+1], y1[2*i+3], yd1[2*i+1], yd1[2*i+3])
            # 2. apply it to the lowest mass (the mass next to the kite)   
            spring_forces = force - last_force    
            last_force = force       
            mass = calc_particle_mass(i+2)
            acc1 = spring_forces / mass  # create the vector of the spring acceleration    
            res_4 = yd1[2*i+4] - (G_EARTH - acc1) # the derivative of the velocity must be equal to the total acceleration  
            RESULT[2*i+3] = res_3 
            RESULT[2*i+4] = res_4                       
    
        # 3. calculate the force of the spring above    
        res_1   = y1[2]  - yd1[1] # the derivative of the position of mass1 must be equal to its velocity
        force = calc_spring_force(y1[0], y1[1], yd1[0], yd1[1])
    
        # 2. apply it to the next mass nearer to the winch
        spring_forces = force - last_force    
        mass = calc_particle_mass(1)
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
    dt = 0.15
    plt.ioff()
    fig = plt.figure()
    plt.ylim(-SEGMENTS*L0-10, 0.5)
    plt.xlim(-SEGMENTS*L0/2, SEGMENTS*L0/2)
    plt.grid(True, color="grey", linestyle="dotted")
    line, sc, txt = None, None, None
    for t in np.linspace(0, duration, num=round(duration/dt)+1):
        line, sc, txt = plot2d(fig, t_sol, y, t, SEGMENTS, line, sc, txt)
        time.sleep(dt/2)
    if os.environ.get("TETHERS_BRIEF_PLOT") == "1":
        # show briefly and close automatically, e.g. when running the tests
        plt.pause(1)
        plt.close('all')
    else:
        plt.show()

def run_example():
    # Create an instance of the problem
    model = ExtendedProblem()  # Create the problem
    model.name = 'Mass-Spring' # Specifies the name of problem (optional)

    sim = IDA(model) # Create the solver
    sim.verbosity = 30
    sim.atol = 1.0e-6
    sim.rtol = 1.0e-6

    t_sol, y, yd = sim.simulate(DURATION, round(DURATION*50)) # 50 communication points per second

    # extract the z position and velocity of the lowest mass (mass SEGMENTS)
    pos_z_ix = 5 + (SEGMENTS - 1) * 6
    vel_z_ix = pos_z_ix + 3
    pos_z = y[:, pos_z_ix]
    vel_z = y[:, vel_z_ix]

    # saving the result for comparison with the Julia implementation
    os.makedirs("output", exist_ok=True)
    with open(os.path.join("output", "Tether_05_python.csv"), "w") as f:
        f.write("time,pos_z,vel_z\n")
        for t_i, pz_i, vz_i in zip(t_sol, pos_z, vel_z):
            f.write(f"{t_i},{pz_i},{vz_i}\n")

    play(DURATION, t_sol, y)
    
if __name__ == '__main__':
    run_example()
