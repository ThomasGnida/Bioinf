"""
Created on Tuesday Mar 24 10:27 2026
Quicksort Assignment IBSY 2026

@author: Gnida
"""
import numpy as np
import random as rnd
import matplotlib.pyplot as plt
from time import time
from math import floor

'''
Code to compare the sorting algorithms was reused as provided by the lecturer.
Small adjustments were made:
    Changed Emin and Emax for myqucksort to mimic tests for built-in quicksort.
    Changed np.sort kind to quicksort 
    Added a legend for the final plot.
    
'''
def check_my_sort():
    # general initialization
    N_pts = 5  # number of different lentghs for which sort time is measured
    Nt = 1  # number of repetitions to measure time
    results = np.zeros([2, N_pts, 2])

    # range of lengths for built in method
    Emin = 5  # log10 of minimal sample size
    Emax = 6  # log10 of maximal sample size
    pts = np.linspace(Emin, Emax, N_pts)  # contains exponents for list lengths
    for i in range(N_pts):
        N = int(10 ** (pts[i]))  # lentgh of list
        num_list = np.zeros(N)
        # inititalization for timing of list with length N
        for j in range(N):
            num_list[j] = rnd.uniform(0, 100)
        # timing for built in method
        t1 = time()
        for t in range(Nt):
            sorted_num_list = np.sort(num_list, kind='quicksort')
            #Changed to quicksort
        t2 = time()
        # write time for each list length in results
        results[0, i, 0] = N
        results[0, i, 1] = (t2 - t1) / Nt

    Emin = 5  # log10 of minimal sample size
    Emax = 6  # log10 of maximal sample size
    pts = np.linspace(Emin, Emax, N_pts)  # contains exponents for list lengths
    for i in range(N_pts):
        N = int(10 ** (pts[i]))  # lentgh of list
        num_list = np.zeros(N)
        # inititalization for timing of list with length N
        for j in range(N):
            num_list[j] = rnd.uniform(0, 100)
        # timing for my clumsy sort
        t1 = time()
        for t in range(Nt):
            sorted_num_list = myquicksort(num_list)
        t2 = time()
        # write time for each list length in results
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

    #Added Legend
    plt.figure()
    plt.plot(x1, y1, 'x' )
    plt.plot(x1, k1 * x1 + d1, '-', label='Numpy Quicksort ')
    plt.ylabel('log(t)')
    plt.xlabel('log(N)')
    print('Slope and intercept of numpy quicksort in log-log plot are ' + str(k1)
          + ' and ' + str(d1))

    plt.plot(x2, y2, 'x')
    plt.plot(x2, k2 * x2 + d2, '-', label='Myquicksort')
    print('Slope and intercept of myquicksort in log-log plot are ' + str(k2)
          + ' and ' + str(d2))
    plt.legend(loc='upper left')
    plt.show()

def myquicksort(list):
    """
    Custom quick sort algorithm
    :param list: List of numbers
    :return: Sorted list of numbers
    """
    if len(list) <= 1:
        return list
        #If there is only one or no element in the array
        # it is sorted and can be returned

    # Median array index is selected as a pivot index
    pivot_index = (floor(len(list)/2))
    pivot = list[pivot_index]
    left = []
    right = []
    for index in range(len(list)):
        # Pivot index is ignored as it will be sorted between the two arrays
        if index == pivot_index:
            continue
        if list[index] < pivot:
            left.append(list[index])
        else:
            right.append(list[index])

    # Function is called recursively
    return myquicksort(left) + [pivot] +  myquicksort(right)



example_array = [0.1, 0.2, 15.2, 10.2, 1.0 , 1.4, 9.3, 13.2, 100.2, 151231.01]
print("Original list:")
print(example_array)
example_array = myquicksort(example_array)
print("Sorted list:")
print(example_array)
check_my_sort()


