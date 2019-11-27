"""
Dual-Population Genetic Algorithm (dgpa) for Adaptive Diversity Control
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
    """
    return 1.00001*max(M) - rosenbrock(M) 

def fr(R):
    """
    Same for the reserve population.
    """
    frvalues = []
    for i in range(R):
        frvalues.append(
    return 1 - abs(10-

def dpga(tmax):
    """
    Maintain population diversity with genetic algorithms. Dual-population
    genetic algorithm does that using an extra population for more diversity.
    Loop until tmax.
    """
    for i in range(tmax):
        GG  
