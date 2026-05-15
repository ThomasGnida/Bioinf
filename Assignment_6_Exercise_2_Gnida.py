"""
Created on Friday May 15 15:07 2026
ODE Assignment IBSY 2026

@author: Gnida
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.special.cython_special import k1


#RK4_step and RK4integrator are reused from exercise 1.
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

    y0 = np.array(y0)

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


def f(t,y):
    #Params derived from fig2 and fig3
    c,m,x = y
    vi = 0.025
    vd = 0.25
    Kd = 0.02
    kd = 0.01
    kc = 0.5
    k = 0.005
    k1 = k2 = k3 = k4 = k

    #Params from fig2 and equation [2]
    vm1 = 3
    vm3 = 1
    v1 = (c / (kc + c)) * vm1
    v2 = 1.5
    v3 = m * vm3
    v4 = 0.5

    #Equations from equation [1]
    dc = vi - vd * x * (c / (Kd + c)) - kd * c
    dm = v1 * (1 - m) / (k1 + (1 - m)) - v2 * m / (k2 + m)
    dx = v3 * (1 - x) / (k3 + (1 - x)) - v4 * x / (k4 + x)

    return np.array([dc,dm,dx])

def exercise_2():
    #Initialization
    t0 = 0
    t_end = 100
    dt = 0.01
    N = int((t_end - t0) / dt)
    y0 = [0.01,0.01,0.01]

    #Calling the rk4 integrator with the previously defined equations
    t, y = RK4integrator(y0, t0, N, f, dt)
    c = y[:, 0]
    m = y[:, 1]
    x = y[:, 2]

    #plotting the figure as close as possible to the original
    fig, ax1 = plt.subplots()

    ax2 = ax1.twinx()
    #Original black and white was changed for improved readability
    ax1.plot(t, m, label="cdc2 kinase (M)", linewidth=2)
    ax1.plot(t, x,'--', label="Protease (X)", linewidth=2, )
    ax2.plot(t, c, label="Cyclin (C)", linewidth=2, color = "green")

    plt.text(16, 0.48, "C", fontsize=14, fontweight='bold')
    plt.text(23, 0.68, "M", fontsize=14, fontweight='bold')
    plt.text(22.5, 0.08, "X", fontsize=14, fontweight='bold')

    ax1.set_ybound(0,1)
    ax1.set_xbound(0,100)
    ax1.xaxis.set_minor_locator(ticker.AutoMinorLocator(4))
    ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax1.tick_params(which='minor', length=3)
    ax1.tick_params(which='major', length=6)
    ax1.tick_params(which='both', direction='in', top=True, right=True)
    ax2.set_ybound(0,1)
    ax2.tick_params(which='major', length=6, direction='in')

    plt.title("Exercise 2 \n Mitotic Oscillator Cascade Model")
    ax1.set_xlabel('Time (min)',fontweight = 'bold')
    ax1.set_ylabel('Fraction of active cdc2 kinase (M)\n'
                   'or cyclin protease (X)',fontweight = 'bold')
    ax2.set_ylabel('Cyclin concentration, C(µM)', fontweight = "bold")
    plt.show()

    return

exercise_2()