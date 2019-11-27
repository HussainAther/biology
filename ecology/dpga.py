import numpy as np

"""
Dual-Population Genetic Algorithm (dgpa) 
Park, Taejin and and Kwang Ryel Ryu "A Dual-Population Genetic Algorithm
for Adaptive Diversity Control" IEEE Transactions on Evolutionary Computation,
Vol 14. No. 6, Dec 2010.
"""

m = 10 # main population initial amount
n = 3 # reserve population initial amount

M = [m] # list of main population amounts at each time interval
R = [n] # same for reserve pop

def rosenbrock(X):
    """
    This R^2 -> R^1 function should be compatible with algopy.
    http://en.wikipedia.org/wiki/Rosenbrock_function
    A generalized implementation is available
    as the scipy.optimize.rosen function.
    """
    x = X[0]
    y = X[1]
    a = 1. - x
    b = y - x*x
    return a*a + b*b*100.

def fm(M):
    """
    Use a fitness function to evaluate the main population.
    The paper tested many functions. Here, we use the Rosenbrock
    function to show the control.
    """
    return 1.00001*max(M) - rosenbrock(M) 

def fr(R, M, sigma=.5):
    """
    Same for the reserve population for sigma distance between the
    two populations.
    """
    frvalues = [] # fitness values
    Mavg = sum(M)/len(M) # average of the main population
    for i in range(R):
        # Calculate each fitness value with the sigma and Euclidean
        # distance. The paper goes into more distances such as the normalized
        # average binary distance.
        frvalues.append(1 - abs(sigma - np.linalg.norm(M-i)))
    return frvalues

def dpga(tmax):
    """
    Maintain population diversity with genetic algorithms. Dual-population
    genetic algorithm does that using an extra population for more diversity.
    Loop until tmax.
    """
    for i in range(tmax):
          
