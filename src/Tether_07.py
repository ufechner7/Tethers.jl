# -*- coding: utf-8 -*-
"""
Tutorial example simulating a 3D mass-spring system with a nonlinear spring (1%
stiffness for l < l_0), five tether segments, tether drag and reel-out.
The DAE is solved with Assimulo's IDA. This is the Python counterpart of
Tether_07.jl, extending Tether_06.py with an aerodynamic drag force per segment.
"""
import numpy as np
import matplotlib.pyplot as plt
import math
import os
import time

from assimulo.solvers.sundials import IDA     # Imports the solver IDA from Assimulo
from assimulo.problem import Implicit_Problem # Imports the problem formulation from Assimulo

G_EARTH       = np.array([0.0, 0.0, -9.81])   # gravitational acceleration        [m/s²]
V_WIND_TETHER = np.array([2.0, 0.0, 0.0])     # wind velocity acting on the tether  [m/s]
RHO           = 1.225                          # air density                     [kg/m³]
CD_TETHER     = 0.958                          # drag coefficient of the tether       [-]
C_SPRING      =  614600.0                      # spring constant
DAMPING       =  473                           # unit damping constant [Ns], must match Tether_07.jl
REL_COMPRESSION_STIFFNESS = 0.01               # relative compression stiffness       [-]
L0            =  50.0                          # initial segment length     [m]
ALPHA0        = math.pi/10                     # initial tether angle     [rad]
SEGMENTS      = 5
DURATION      = 10                             # duration of the simulation [s]

V_RO = 2.0                             # reel-out speed                  [m/s]
D_TETHER = 4.0                         # tether diameter                  [mm]
RHO_TETHER = 724.0                     # density of Dyneema            [kg/m³]

ZEROS  = np.array([0.0, 0.0, 0.0])
RESULT = np.zeros(SEGMENTS * 6 + 3).reshape((-1, 3))
mass_per_meter = RHO_TETHER * math.pi * (D_TETHER/2000.0)**2

def calc_c_spring(norm, l_seg, c_spring0):
    """ Nonlinear spring stiffness of a segment: full stiffness while taut
        (norm > l_seg), only REL_COMPRESSION_STIFFNESS of it while slack. Unlike
        Tether_06.py, the switch is a hard step, matching Tether_07.jl. """
    taut = 1.0 if norm > l_seg else 0.0
    return c_spring0 / (1.0 + REL_COMPRESSION_STIFFNESS) * (REL_COMPRESSION_STIFFNESS + taut)

def calc_spring_force(pos1, pos2, vel1, vel2, l_seg, c_spring0, damping):
    """ Spring and damping force of the segment between the masses at pos1 and pos2.
        The result points from pos1 towards pos2. Spring and damper act in parallel
        along the segment, therefore the damping force uses the component of the
        relative velocity along the segment and not the full 3D vector. """
    segment     = pos2 - pos1
    norm        = np.linalg.norm(segment)
    unit_vector = segment / norm
    c_spring    = calc_c_spring(norm, l_seg, c_spring0)
    spring_vel  = np.dot(vel2 - vel1, unit_vector)   # rate of change of the segment length
    return (c_spring * (norm - l_seg) + damping * spring_vel) * unit_vector

def calc_drag_force(pos1, pos2, vel1, vel2):
    """ Half of the aerodynamic drag force of the segment between pos1 and pos2,
        matching half_drag_force in Tether_07.jl. Unlike the spring force, this is
        not an action/reaction pair: the same vector is added to both end masses,
        since it is one half of the drag of the segment they share. """
    segment     = pos2 - pos1
    norm        = np.linalg.norm(segment)
    unit_vector = segment / norm
    v_apparent  = V_WIND_TETHER - (vel1 + vel2) / 2.0
    v_app_perp  = v_apparent - np.dot(v_apparent, unit_vector) * unit_vector
    norm_v_app  = np.linalg.norm(v_app_perp)
    return 0.25 * RHO * CD_TETHER * norm_v_app * (norm * D_TETHER / 1000.0) * v_app_perp

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
    # Set the initial conditions. Unlike Tether_06.jl, the initial velocities already
    # follow the reel-out speed (consistent with a tether growing at V_RO), matching
    # calc_initial_state in Tether_07.jl, which reduces the initial transient.
    t0  = 0.0                   # Initial time
    pos, vel, acc = [], [], []
    for i in range (SEGMENTS + 1):
        l0_ = -i*L0/SEGMENTS
        v0_ = -i*V_RO/SEGMENTS
        pos.append(np.array([math.sin(ALPHA0) * l0_, 0.0, math.cos(ALPHA0) * l0_]))
        vel.append(np.array([math.sin(ALPHA0) * v0_, 0.0, math.cos(ALPHA0) * v0_]))
        if i == 0:
            acc.append(np.array([0.0, 0.0, 0.0]))
        else:
            acc.append(G_EARTH)
    y0, yd0 = pos[0], vel[0]
    for i in range (SEGMENTS):
        y0  = np.append(y0,  np.append(pos[i+1], vel[i+1])) # Initial state vector
        yd0 = np.append(yd0, np.append(vel[i+1], acc[i+1])) # Initial state vector derivative
    # mass0's position is a pure algebraic constraint (res only involves y, never yd
    # for these three components), the rest are true differential states
    algvar = np.append(np.zeros(3), np.ones(6 * SEGMENTS))

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
        last_drag  = ZEROS

        for i in range(SEGMENTS-2, -1, -1):    # count down from segments-2 to zero
            # 1. calculate the force of the lowest spring (the spring next to the kite)
            res_3   =  y1[2*i+4] - yd1[2*i+3]  # the derivative of the position of mass1 must be equal to its velocity
            # the spring and drag force of the segment between the masses i+1 and i+2
            p1, p2, v1, v2 = y1[2*i+1], y1[2*i+3], yd1[2*i+1], yd1[2*i+3]
            force = calc_spring_force(p1, p2, v1, v2, l_seg, c_spring1, damping)
            drag  = calc_drag_force(p1, p2, v1, v2)
            # 2. apply them to the lowest mass (the mass next to the kite). The spring
            # force is an action/reaction pair, so the two segments of a mass cancel
            # each other, while both segments push it downwind with their drag half.
            spring_forces = force - last_force
            drag_forces   = drag + last_drag
            last_force = force
            last_drag  = drag
            mass = calc_particle_mass(i+2, m_tether_particle)
            # acc1 enters the residual below as (G_EARTH - acc1), therefore the drag,
            # which acts along the apparent wind, has to be subtracted here
            acc1 = (spring_forces - drag_forces) / mass
            res_4 = yd1[2*i+4] - (G_EARTH - acc1) # the derivative of the velocity must be equal to the total acceleration
            RESULT[2*i+3] = res_3
            RESULT[2*i+4] = res_4

        # 3. calculate the force of the spring above
        res_1   = y1[2]  - yd1[1] # the derivative of the position of mass1 must be equal to its velocity
        force = calc_spring_force(y1[0], y1[1], yd1[0], yd1[1], l_seg, c_spring1, damping)
        drag  = calc_drag_force(y1[0], y1[1], yd1[0], yd1[1])

        # 2. apply it to the next mass nearer to the winch
        spring_forces = force - last_force
        drag_forces   = drag + last_drag
        mass = calc_particle_mass(1, m_tether_particle)
        acc = (spring_forces - drag_forces) / mass
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
    if os.environ.get("TETHERS_BRIEF_PLOT") == "1":
        # show briefly and close automatically, e.g. when running the tests
        plt.pause(1)
        plt.close('all')
    else:
        plt.show(block=True)

def run_example():
    # Create an instance of the problem
    model = ExtendedProblem()  # Create the problem
    model.name = 'Mass-Spring' # Specifies the name of problem (optional)

    sim = IDA(model) # Create the solver, using the default dense direct linear solver
    sim.verbosity = 0
    sim.atol = 1.0e-5
    sim.rtol = 1.0e-5
    sim.algvar = model.algvar
    sim.suppress_alg = True
    sim.maxord = 3

    start = time.perf_counter()
    t_raw, y_raw, _ = sim.simulate(DURATION, round(DURATION*50)) # 50 communication points per second
    elapsed_time = time.perf_counter() - start
    print(f"Elapsed time: {elapsed_time} s, speed: {round(DURATION/elapsed_time)} times real-time")

    # IDA's chosen output points can silently land off the shared 0.02s grid. Resample the
    # full state onto the exact grid Tether_07.jl uses (ts=0:dt:duration), so the two CSVs
    # are directly comparable and play() gets a consistent (t_sol, y) pair.
    t_sol = np.linspace(0.0, DURATION, round(DURATION*50) + 1)
    y = np.column_stack([np.interp(t_sol, t_raw, y_raw[:, j]) for j in range(y_raw.shape[1])])

    # extract the z position and velocity of the lowest mass (mass SEGMENTS)
    pos_z_ix = 5 + (SEGMENTS - 1) * 6
    vel_z_ix = pos_z_ix + 3
    pos_z = y[:, pos_z_ix]
    vel_z = y[:, vel_z_ix]

    # saving the result for comparison with the Julia implementation
    os.makedirs("output", exist_ok=True)
    with open(os.path.join("output", "Tether_07_python.csv"), "w") as f:
        f.write("time,pos_z,vel_z\n")
        for t_i, pz_i, vz_i in zip(t_sol, pos_z, vel_z):
            f.write(f"{t_i},{pz_i},{vz_i}\n")

    play(DURATION, t_sol, y)
    return

if __name__ == '__main__':
    run_example()
