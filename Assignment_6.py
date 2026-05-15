"""
. IVP Solver
Write a single python module that:
A) can solve the general initial value problem (IVP):
𝒚′ = 𝒇(𝑡, 𝒚) ; 𝒚(𝑡0
) = 𝒚0
using the explicit 4-th order Runge-Kutta method (RK4). (Note: y, y0 and f are vectors!)
B) solves the following concrete IVP in the interval [0,10] with your RK4 implementation:
𝑦
′′ = −4𝑦 − 𝑦
′
𝑦(0) = 1; 𝑦
′
(0) = −1
and compares the result with the RK45 solution obtained from
scipy.integrate.solve_ivp.1

The module should contain the following:
• A function RK4_step(t,y,dt,f) that performs a single RK4 integration step and
returns the next time and (approximate) function values (t_next,y_next).
𝑅𝐾4_𝑠𝑡𝑒𝑝 (𝑡𝑛, 𝑦𝑛, 𝑑𝑡, 𝑓) → 𝑟𝑒𝑡𝑢𝑟𝑛 (𝑡𝑛+1, 𝑦𝑛+1)
The input arguments are:
o t: current time (scalar)
o y: current function values (vector)
o dt: time step (scalar)
o f: function to calculate the derivative (function).
Note: it should work for an arbitrary function
𝑓: ℝ × ℝ𝑛 → ℝ𝑛
; (𝑡, 𝒚) → 𝒇(𝑡, 𝒚)


• A function RK4integrator(y0,f,N,t0,dt) that performs:
o N integration steps
1 RK45 (Dormand-Prince [2], [3]) is the default method of scipy.integrate.solve_ivp.
https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html
3
o of length dt
o starting at initial conditions (t0,y0)
by calling the RK4_step function N-times. It should return (t,Y):
o t: an array containing the time grid
o Y: an array containing the numerical solution to the general IVP.
𝑌[𝑛, :] = 𝒚𝑛


• A function test_integration() for task B. It should plot both solutions in a single
figure.


2. Minimal Cascade Model
Write a python module that reproduces figure 3 of reference [1]. The output could look like
this:


"""

from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

def RK4_step(t, y, dt, f):
    k1 = f(t,y)
    k2 = f(t + dt/2, y + dt/2 * k1)
    k3 = f(t + dt/2, y + dt/2 * k2)
    k4 = f(t + dt, y + dt * k3)
    y_next = y + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
    return (t + dt, y_next)

def RK4_integrator(y0, f, N, t0, dt):
    y0 = np.array(y0, dtype = float)
    t = np.zeros(N+1)
    y = np.zeros((N+1, len(y0)))

    t[0] = t0
    y[0, :] = y0

    curent_t = t0
    curent_y = y0

    for n in range (N):
        curent_t, curent_y = RK4_step(curent_t, curent_y, dt, f)
        t[n+1] = curent_t
        y[n+1, :] = curent_y

    return t, y

test_integration()


