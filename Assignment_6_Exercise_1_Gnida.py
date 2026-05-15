"""
Created on Thursday May 14 15:04 2026
ODE Assignment IBSY 2026

@author: Gnida
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


def RK4_step(t, y, dt, f):
    """
    Perform a single Runge Kutter 4 (RK4) integration step

    :param  t: Current time
    :param  y: Current function values
    :param dt: Time step
    :param  f: Derivative function f(t, y)
    :return: t_next, y_next = next timestamp and solution values
    """

    k1 = f(t, y)
    k2 = f(t + dt/2, y + dt*k1/2)
    k3 = f(t + dt/2, y + dt*k2/2)
    k4 = f(t + dt, y + dt*k3)

    y_next = y + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)
    t_next = t + dt

    return t_next, y_next


def RK4integrator(y0, t0 , N, f, dt):
    """
    Perform N RK4 integration steps.

    :param y0: Initial condition
    :param t0: Initial time
    :param N: Number of steps
    :param f: Derivative function f(t, y)
    :param dt: Time step
    :return: t, Y = time grid and numerical solution
    """

    y0 = np.array(y0, dtype=float)

    t = np.zeros(N + 1)
    y = np.zeros((N + 1, len(y0)))

    t[0] = t0
    y[0] = y0

    current_t = t0
    current_y = y0

    for n in range(N):
        current_t, current_y = RK4_step(current_t, current_y, dt, f)

        t[n + 1] = current_t
        y[n + 1] = current_y

    return t, y

def test_integration():
    '''
    Uses the predefined parameters to compare rk4 and rk45 from scipy.
    No return but prints the results for t_end and plots both methods in a
    matplotlib figure.
    '''

    t_start = 0
    t_end = 10
    y0 = np.array([1.0, -1.0])
    dt = 0.01
    N = int((t_end - t_start) / dt)

    #Define function as first order system
    def f(t, y):
        y1, y2 = y
        return np.array([y2, -4*y1 - y2 ])

    t_rk4, y_rk4 = RK4integrator(y0, t_start, N, f, dt)
    #Call the solve_ivp function with t_rk4 as evaluation timesteps
    rk45_solution = solve_ivp(f,
                              (t_start, t_end),
                              y0,
                              t_eval=t_rk4,
                              method='RK45') #Default = RK45

    #Figure creation
    plt.figure(figsize=(10, 6))

    plt.plot(t_rk4, y_rk4[:, 0], label='RK4: y(t)', linewidth=3)
    plt.plot(rk45_solution.t,
             rk45_solution.y[0],
             '--',
             label='RK45: y(t)',
             linewidth=2)
    plt.plot(t_rk4, y_rk4[:, 1], label="RK4: y'(t)", linewidth=3)
    plt.plot(rk45_solution.t,
             rk45_solution.y[1],
             '--',
             label="RK45: y'(t)",
             linewidth=2)

    plt.title('Task B: Comparison of RK4 and RK45 Solutions')
    plt.xlabel('t')
    plt.ylabel('Approximate Solution')
    plt.grid(True)
    plt.legend()
    plt.show()

    #Printing final values for comparison
    print("Final RK4 solution:")
    print("y(10) =", y_rk4[-1, 0])
    print("y'(10)=", y_rk4[-1, 1])

    print("\nFinal RK45 solution:")
    print("y(10) =", rk45_solution.y[0, -1])
    print("y'(10)=", rk45_solution.y[1, -1])
    return

test_integration()