# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 18:45:10 2020
Demonstration of scaling behavior with sort algorithms, using a log-log plots.

@author: Ramberger
"""

import numpy as np
import random as rnd
import matplotlib.pyplot as plt
from time import time

# a personal (clumsy) sort alogrithm that iterates the current position through
# the whole array (outer loop i), and then always swaps a smaller element with 
# the current element if it finds one (inner loop j)
def my_sort(list_to_sort):
    N = len(list_to_sort)
    list_sorted = list_to_sort.copy()
    
    for i in range(N):
        for j in range(i,N-1): #would be even slower, if we always started at 0
            if list_sorted[i] > list_sorted[j + 1]:
                list_sorted[i] , list_sorted[j + 1] = \
                list_sorted[j + 1], list_sorted[i]
    return list_sorted

#Compare scaling of different sort algos
def check_my_sort():
    #general initialization
    N_pts = 5 # number of different lentghs for which sort time is measured
    Nt = 1 #number of repetitions to measure time
    results = np.zeros([2, N_pts, 2])
    
    #range of lengths for built in method
    Emin = 5  # log10 of minimal sample size
    Emax = 6  # log10 of maximal sample size
    pts = np.linspace(Emin, Emax, N_pts) #contains exponents for list lengths
    for i in range(N_pts):
        N = int(10 ** (pts[i])) #lentgh of list
        num_list = np.zeros(N)
        print(N)
        #inititalization for timing of list with length N
        for j in range(N):
            num_list[j] = rnd.uniform(0, 100)
        #timing for built in method
        t1 = time()
        for t in range(Nt):
            sorted_num_list = np.sort(num_list, kind='mergesort')
        t2 = time()
        #write time for each list length in results
        results[0, i, 0] = N
        results[0, i, 1] = (t2 - t1) / Nt
        
    #range of lengths for my clumsy sort
    Emin = 2  # log10 of minimal sample size
    Emax = 3  # log10 of maximal sample size
    pts = np.linspace(Emin, Emax, N_pts) #contains exponents for list lengths
    for i in range(N_pts):
        N = int(10 ** (pts[i])) #lentgh of list
        num_list = np.zeros(N)
        print(N)
        #inititalization for timing of list with length N
        for j in range(N):
            num_list[j] = rnd.uniform(0, 100)
        #timing for my clumsy sort
        t1 = time()
        for t in range(Nt):
            sorted_num_list = my_sort(num_list)
        t2 = time()
        #write time for each list length in results
        results[1, i, 0] = N
        results[1, i, 1] = (t2 - t1) / Nt

    print(results[0, :, :])
    print(results[1, :, :])
    x1 = np.log10(results[0, :, 0])
    y1 = np.log10(results[0, :, 1])
    x2 = np.log10(results[1, :, 0])
    y2 = np.log10(results[1, :, 1])
    k1, d1 = np.polyfit(x1, y1, 1)
    k2, d2 = np.polyfit(x2, y2, 1)

    plt.figure()
    plt.plot(x1, y1, 'x')
    plt.plot(x1, k1 * x1 + d1, '-')
    plt.ylabel('log(t)')
    plt.xlabel('log(N)')
    print('Slope and intercept of built-in in log-log plot are ' + str(k1)
          + ' and ' + str(d1))

    plt.plot(x2, y2, 'x')
    plt.plot(x2, k2 * x2 + d2, '-')
    print('Slope and intercept of my clumsy sort in log-log plot are ' + str(k2)
          + ' and ' + str(d2))
    plt.show()

if __name__ == "__main__":
    check_my_sort()