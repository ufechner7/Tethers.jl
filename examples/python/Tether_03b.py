# -*- coding: utf-8 -*-
"""
Tutorial example simulating a falling mass, attached to a non-linear spring
(no spring force while the segment is loose), with a state event at the transition.

The model is written as a CasADi expression graph and integrated with SUNDIALS' CVODES,
which CasADi supplies with the exact Jacobian of the model.

The taut/slack crossings are located by hand, which is where Python is genuinely behind
Julia here. CasADi does have event detection - a `zero` entry in the DAE dictionary - but
it is experimental, and on this model it aborts part way through with "tout too far back
in direction of integration": after restarting at an event it can no longer interpolate
back to the requested output point. So `run_example` walks the output grid itself,
watches the sign of the segment length minus its rest length, bisects for the crossing
and restarts the integrator there. No state has to change at the crossing: the switch is
recomputed from the current position on every call. Placing an exact step boundary there
is the whole point, and it is what `Tether_03b.jl` gets from a `ContinuousCallback` and
what Assimulo gets from `state_events`.

To advance by an arbitrary amount without rebuilding the integrator, the time is scaled:
the solver integrates `tau` from 0 to 1 over an ODE multiplied by the step length `h`,
which is passed as a parameter.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import casadi as ca

G_EARTH  = np.array([0.0, 0.0, -9.81]) # gravitational acceleration
C_SPRING =  50                         # spring constant [N/m]; Dynema, 4mm: 200e3
DAMPING  =  0.5                        # damping [Ns/m]
MASS     = 1.0                         # mass per point-mass [kg]
L_0      = 10.0                        # initial segment length [m]
V0       = 4.0                         # initial velocity

# Falling mass, attached to a spring anchored at the origin
# State vector y = mass0.pos, mass1.pos, mass1.vel
def build_model():
    """ The two masses and their spring as one CasADi expression graph. Mass 0 is fixed,
        so its derivative is zero and it stays at the origin. Returns the state vector
        `y`, its derivative `ydot`, the event indicator and the initial state `y0`. """
    y = ca.SX.sym('y', 9)
    segment = y[3:6] - y[0:3]       # the vector from mass0 to mass1
    norm = ca.norm_2(segment)
    unit_vector = segment / norm
    rel_vel = y[6:9]                # mass0 is fixed, so this is mass1's velocity
    # if the segment is loose (norm <= L_0) there is no spring force at all
    c_spring = ca.if_else(norm > L_0, C_SPRING, 0.0)
    # spring and damper act in parallel along the segment, therefore the damping
    # force uses the component of the relative velocity along the segment
    spring_vel = ca.dot(rel_vel, unit_vector)
    force = (c_spring * (norm - L_0) + DAMPING * spring_vel) * unit_vector
    acc = force / MASS
    ydot = ca.vertcat(ca.SX.zeros(3), y[6:9], G_EARTH - acc)
    y0 = np.array([0.0, 0.0, 0.0, 0.0, 0.0, -L_0, 0.0, 0.0, V0])
    return y, ydot, norm - L_0, y0


def c_spring_of(y_sol):
    """ Spring constant at every sample, for the grey line in the plot. """
    norm = np.linalg.norm(y_sol[:, 3:6] - y_sol[:, 0:3], axis=1)
    return np.where(norm > L_0, C_SPRING, 0.0)


def build_integrator():
    """ An integrator that advances the state by an arbitrary step length `h`, and a
        function evaluating the event indicator. Time is scaled so that the solver always
        integrates tau from 0 to 1, with the step length as a parameter. """
    y, ydot, event, y0 = build_model()
    h = ca.SX.sym('h')
    sim = ca.integrator('sim', 'cvodes', {'x': y, 'p': h, 'ode': h * ydot}, 0.0, 1.0,
                        {'abstol': 1.0e-6, 'reltol': 1.0e-6})
    return sim, ca.Function('indicator', [y], [event]), y0


def run_example():
    sim, indicator, y0 = build_integrator()

    def step(state, dt):
        return np.array(sim(x0=state, p=dt)['xf']).flatten()

    time = np.linspace(0.0, 10.0, 501)
    states, current, crossings = [y0], y0, 0
    for k in range(len(time) - 1):
        dt = time[k+1] - time[k]
        before = float(indicator(current))
        after = step(current, dt)
        if before * float(indicator(after)) < 0.0:
            # the segment changed between taut and slack inside this interval: bisect for
            # the crossing and restart the integration exactly there
            lo, hi = 0.0, dt
            for _ in range(40):
                mid = 0.5 * (lo + hi)
                if before * float(indicator(step(current, mid))) < 0.0:
                    hi = mid
                else:
                    lo = mid
            current = step(step(current, hi), dt - hi)
            crossings += 1
        else:
            current = after
        states.append(current)
    print("Number of state events:", crossings)
    y_sol = np.array(states)

    # plot the result
    pos_z = y_sol[:, 5]
    vel_z = y_sol[:, 8]
    C_SPRINGS = c_spring_of(y_sol)

    # saving the result for comparison with the Julia implementation
    os.makedirs("output", exist_ok=True)
    with open(os.path.join("output", "Tether_03b_python.csv"), "w") as f:
        f.write("time,pos_z,vel_z\n")
        for t_i, pz_i, vz_i in zip(time, pos_z, vel_z):
            f.write(f"{t_i},{pz_i},{vz_i}\n")

    plt.gcf().canvas.manager.set_window_title("falling mass, non-linear spring")
    plt.ax1 = plt.subplot(111)
    plt.ax1.set_xlabel('time [s]')
    plt.plot(time, pos_z, color="green")
    plt.plot(time, -np.ones(len(time)) * L_0 + 0.005 * C_SPRINGS, color="grey", label="c_spring")
    plt.ax1.set_ylabel('pos_z [m]')
    plt.grid(True)
    plt.ax2 = plt.twinx()
    plt.ax2.set_ylabel('vel_z [m/s]')
    plt.plot(time, vel_z, color="red")

    if os.environ.get("TETHERS_BRIEF_PLOT") == "1":
        # show briefly and close automatically, e.g. when running the tests
        plt.pause(1)
        plt.close('all')
    else:
        plt.show()


if __name__ == '__main__':
    run_example()
