# -*- coding: utf-8 -*-
"""
Tutorial example simulating a 3D mass-spring system with a nonlinear spring (1% spring
forces for l < l_0, smoothly blended in near l_0), five tether segments and reel-out.
The DAE is solved with Assimulo's IDA, using an analytic Jacobian.
"""
import numpy as np
import matplotlib.pyplot as plt
import math
import os
import time

from assimulo.solvers.sundials import IDA     # Imports the solver IDA from Assimulo
from assimulo.problem import Implicit_Problem # Imports the problem formulation from Assimulo

G_EARTH  = np.array([0.0, 0.0, -9.81]) # gravitational acceleration
C_SPRING =  614600.0                   # spring constant
DAMPING  =  473                        # unit damping constant [Ns], must match Tether_06.jl
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

SMOOTH_REL_WIDTH = 1e-3  # width of the taut/slack blend, as a fraction of l_seg

def _smoothstep(x):
    """ 0 for x<=0, 1 for x>=1, cubic (C1) ramp in between. """
    xc = min(max(x, 0.0), 1.0)
    return xc * xc * (3.0 - 2.0 * xc)

def _smoothstep_deriv(x):
    """ d(_smoothstep)/dx. """
    if x <= 0.0 or x >= 1.0:
        return 0.0
    return 6.0 * x * (1.0 - x)

def calc_c_spring(norm, l_seg, c_spring0):
    """ Spring stiffness of a segment. A taut segment (norm > l_seg) uses the full
        stiffness; a slack one still keeps 1% of it, matching Tether_06.jl. The
        switch between the two is blended smoothly (over SMOOTH_REL_WIDTH of l_seg)
        instead of a hard step, so the force stays C1 and IDA's Newton iteration
        stays well-conditioned right at the transition. """
    width = SMOOTH_REL_WIDTH * l_seg
    frac = _smoothstep((norm - l_seg) / width)
    return c_spring0 / 1.01 * (0.01 + frac)

def calc_c_spring_and_deriv(norm, l_seg, c_spring0):
    """ Same as calc_c_spring, plus d(c_spring)/d(norm). """
    width = SMOOTH_REL_WIDTH * l_seg
    x = (norm - l_seg) / width
    c_spring = c_spring0 / 1.01 * (0.01 + _smoothstep(x))
    dc_dnorm = c_spring0 / 1.01 * _smoothstep_deriv(x) / width
    return c_spring, dc_dnorm

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

def calc_particle_mass(i, m_tether_particle):
    """ Mass of particle i. The last particle is attached to one segment only and
        therefore carries half of the mass of an inner particle. """
    if i == SEGMENTS:
        return 0.5 * m_tether_particle
    return m_tether_particle

def calc_spring_force_jac(pos1, pos2, vel1, vel2, l_seg, c_spring0, damping):
    """ Force of calc_spring_force plus its analytic Jacobians w.r.t. pos1, pos2,
        vel1 and vel2 (each a 3x3 matrix, dF_i/dx_j). """
    segment     = pos2 - pos1
    norm        = np.linalg.norm(segment)
    unit_vector = segment / norm
    c_spring, dc_dnorm = calc_c_spring_and_deriv(norm, l_seg, c_spring0)
    rel_vel     = vel2 - vel1
    spring_vel  = np.dot(rel_vel, unit_vector)
    force_mag   = c_spring * (norm - l_seg) + damping * spring_vel
    force       = force_mag * unit_vector

    # M = d(unit_vector)/d(pos2) = -d(unit_vector)/d(pos1)
    uu = np.outer(unit_vector, unit_vector)
    M  = (np.eye(3) - uu) / norm
    m_relvel = M @ rel_vel  # = d(spring_vel)/d(pos2) = -d(spring_vel)/d(pos1)

    # inside the blend band c_spring itself varies with norm, so its contribution
    # to d(force_mag)/d(norm) is c_spring + dc_dnorm*(norm - l_seg), not just c_spring
    c_eff = c_spring + dc_dnorm * (norm - l_seg)
    d_force_mag_dpos1 = -c_eff * unit_vector - damping * m_relvel
    d_force_mag_dpos2 =  c_eff * unit_vector + damping * m_relvel

    dF_dpos1 = np.outer(unit_vector, d_force_mag_dpos1) - force_mag * M
    dF_dpos2 = np.outer(unit_vector, d_force_mag_dpos2) + force_mag * M
    dF_dvel1 = -damping * uu
    dF_dvel2 =  damping * uu
    return force, dF_dpos1, dF_dpos2, dF_dvel1, dF_dvel2

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

    def jac(self, c, t, y, yd):
        """ Analytic Jacobian J = d(res)/dy + c * d(res)/dyd, replacing the
            finite-difference approximation IDA would otherwise use. """
        y1  = y.reshape((-1, 3))
        yd1 = yd.reshape((-1, 3))
        l_seg = (L0 + V_RO*t) / SEGMENTS
        c_spring1 = C_SPRING / l_seg
        damping  = DAMPING / l_seg
        m_tether_particle = mass_per_meter * l_seg
        n = SEGMENTS
        size = 3 * (2*n + 1)
        J = np.zeros((size, size))
        eye3 = np.eye(3)

        def blk(i):
            return slice(3*i, 3*i + 3)

        def pos_block(k):   # index (in the 3-vector blocks) of the position of mass k
            return 0 if k == 0 else 2*k - 1

        def vel_block(k):   # index of the velocity of mass k (k >= 1)
            return 2*k

        # mass0 is fixed: res = y[pos0]
        J[blk(0), blk(0)] = eye3

        for k in range(1, n + 1):
            pb, vb = pos_block(k), vel_block(k)
            # velocity/position consistency: y[vel_k] - yd[pos_k] = 0
            J[blk(pb), blk(vb)] += eye3
            J[blk(pb), blk(pb)] += -c * eye3

            # acceleration equation: yd[vel_k] - (G - acc) = 0
            J[blk(vb), blk(vb)] += c * eye3
            mass = calc_particle_mass(k, m_tether_particle)

            # force of the segment below (between mass k-1 and mass k), added
            p1b = pos_block(k - 1)
            _, dF_dp1, dF_dp2, dF_dv1, dF_dv2 = calc_spring_force_jac(
                y1[p1b], y1[pb], yd1[p1b], yd1[pb], l_seg, c_spring1, damping)
            J[blk(vb), blk(p1b)] += (dF_dp1 + c * dF_dv1) / mass
            J[blk(vb), blk(pb)]  += (dF_dp2 + c * dF_dv2) / mass

            # force of the segment above (between mass k and mass k+1), subtracted
            if k < n:
                p2b = pos_block(k + 1)
                _, dF_dp1, dF_dp2, dF_dv1, dF_dv2 = calc_spring_force_jac(
                    y1[pb], y1[p2b], yd1[pb], yd1[p2b], l_seg, c_spring1, damping)
                J[blk(vb), blk(pb)]  -= (dF_dp1 + c * dF_dv1) / mass
                J[blk(vb), blk(p2b)] -= (dF_dp2 + c * dF_dv2) / mass

        return J

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
    sim.usejac = True

    t_raw, y_raw, _ = sim.simulate(DURATION, round(DURATION*50)) # 50 communication points per second

    # IDA's chosen output points can silently land off the shared 0.02s grid (e.g. it
    # sometimes drops the point at t=9.98 near the very end). Resample the full state
    # onto the exact grid Tether_06.jl uses (ts=0:dt:duration), so the two CSVs are
    # directly comparable and play() gets a consistent (t_sol, y) pair.
    t_sol = np.linspace(0.0, DURATION, round(DURATION*50) + 1)
    y = np.column_stack([np.interp(t_sol, t_raw, y_raw[:, j]) for j in range(y_raw.shape[1])])

    # extract the z position and velocity of the lowest mass (mass SEGMENTS)
    pos_z_ix = 5 + (SEGMENTS - 1) * 6
    vel_z_ix = pos_z_ix + 3
    pos_z = y[:, pos_z_ix]
    vel_z = y[:, vel_z_ix]

    # saving the result for comparison with the Julia implementation
    os.makedirs("output", exist_ok=True)
    with open(os.path.join("output", "Tether_06_python.csv"), "w") as f:
        f.write("time,pos_z,vel_z\n")
        for t_i, pz_i, vz_i in zip(t_sol, pos_z, vel_z):
            f.write(f"{t_i},{pz_i},{vz_i}\n")

    play(DURATION, t_sol, y)
    return
    
if __name__ == '__main__':
    run_example()
