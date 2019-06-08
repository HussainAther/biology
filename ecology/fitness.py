import sys
import numpy as np

"""
Given a pop, return the best individual and the associated value.
"""

def fitness_function(w, x, y, z):
    """
    Fitness_function of four variables with a maximum at (w=-.25, x=3, y=-1, z=2)
    """
    return -2 * (w ** 2) + np.sqrt(abs(w)) - (x ** 2) + (6 * x) - (y ** 2) - (2 * y) - (z ** 2) + (4 * z)

def simple_fitness_function(x, y):
    """
    Simpler fitness_function of two variables with a maximum at (x=1, y=2)
    """
    return - (x**2) + (2 * x) - (y ** 2) + (4 * y)

def find_best(pop):
    """
    For an input array pop of individual fitness levels, find the best. 
    """
    best = None
    val = None
    for individual in pop:
        if len(individual) == 2:
            r = simple_fitness_function(individual[0], individual[1])
            try:
                if r > val:
                    best = individual
                    val = r
                except:  # On the first run, set the result as best
                    best = individual
                    val = r
        elif len(individual) == 4:
            r = fitness_function(individual[0], individual[1], individual[2], individual[3])
            try:
                if r > val:
                    best = individual
                    val = r
                except:  # On the first run, set the result as best
                    best = individual
                    val = r
        else:
            print("error: Wrong number of arguments received")
