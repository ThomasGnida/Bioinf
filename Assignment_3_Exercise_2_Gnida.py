"""
Created on Wednesday Mar 25 11:38 2026
Monte Carlo Pi Assignment IBSY 2026

@author: Gnida
"""
import random
import numpy as np
import matplotlib.pyplot as plt


def mc_pi(N):
    """
    Approximates pi using a Monte Carlo simulation by creating N points
    inside a 2x2 square and comparing with the unit circle.

    :param N: Number of points to simulate
    :return: Approximation of Pi
    """
    inside = 0
    if N == 0:
        print("N cant be 0")
        return 0
    for i in range(N):
        #Randomly generate N coordinates in our 2x2 square
        x = random.uniform(-1, 1)
        y = random.uniform(-1, 1)

        #If our point is inside the unit circle we advance our counter
        if x ** 2 + y ** 2 <= 1:
            inside += 1

    #As our square is 2x2 we assume that p(inside) = pi/4
    # inside/N ~ pi/4 -->
    pi = 4* inside / N
    return pi


def mc_pi_stat(N,M):
    """
    Perform M independent MC simulation for pi using N points each.
    The function mc_pi(N) is called M times, in order to calculate the mean and
    standard deviation.

    :param N: Number of points to simulate
    :param M: Number of simulations performed
    :return: mean and std of performed simulations
    """
    pi_approximations = []
    for i in range(M):
        pi_approximations.append(mc_pi(N))
    mean = np.mean(pi_approximations)
    std = np.std(pi_approximations,ddof=1)
    return mean, std


def mc_pi_plt():
    """
    Performs MC pi estimations for multiple sample sizes and visualizes
    the results.
    Steps:
    1. Uses predefined sample sizes to estimate pi with mc_pi_stat.
    2. Means and std returned from mc_pi_stat are plotted as error bars
        against log10(sample size).
    3. Plots the change of standard deviation std(N) against log10(Sample size).
    4. Calculates and prints the scaling exponent a using std(N) ~ N^a

    :return: No returns,
        the function only produces 2 plots and a print statement
    """
    sample_sizes = [10,32,100,320,1000,10000,100000]
    std_logs = []

    #Plot 1: Predicted Pi for log(N) samples
    plt.title("Monte Carlo Pi Estimation")
    plt.xlabel("log(N)")
    plt.ylabel("µ(N)")
    plt.axhline(np.pi, linestyle='--', label="Pi") #Plot exact pi

    for N in sample_sizes:
        #Call mc_pi_stat for each sample size N simulating 10 coordinates.
        #Log10 of std is stored in an array for future calculations, mean and
        # std are plotted as pyplot.errorbar
        mean, std = mc_pi_stat(N,10)
        std_logs.append(np.log10(std))
        plt.errorbar(np.log10(N), mean, yerr=std, fmt='o', label=f"N={N}")
    plt.legend(ncols=2)
    plt.show()


    #Plot 2: Scaling of the error
    plt.title("Scaling of error")
    plt.xlabel("log(N)")
    plt.ylabel("log(std[N]))")

    #Convert sample-sizes into log10 numpy array and
    # plot log10(std) against log10(N)
    sample_sizes = np.array(np.log10(sample_sizes))
    plt.plot(sample_sizes, std_logs, 'o-')
    plt.show()



    # 3 Compute the scaling coefficient alpha.
    # This is done by calculating the difference between neighbouring points and
    # averaging them.

    alphas = []
    # Compute alpha between each pair of neighboring points
    for i in range(len(sample_sizes) - 1):
        d_std_log = std_logs[i + 1] - std_logs[i]
        d_log_N = sample_sizes[i + 1] - sample_sizes[i]
        alphas.append(d_std_log / d_log_N)

    alpha = np.mean(alphas)
    print("Estimated alpha =", alpha)


mc_pi_plt()
'''
As we expect alpha to be -0.5 the results are confirming our expectations.
To reduce the std by 50% we would need to increase N by a factor of 4.
'''