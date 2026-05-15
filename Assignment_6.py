"""
Created on Thursday May 14 15:04 2026
Quicksort Assignment IBSY 2026

@author: Gnida
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


def RK4_step(t, y, dt, f):
    """
    Perform a single Runge Kutter 4 (RK4) integration step

    :param  t: Current time
    :param  y: Current solution
    :param dt: Time step
    :param  f: Derivative function f(t, y)
    :return: t_next, y_next = next timestamp and solution vector
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
