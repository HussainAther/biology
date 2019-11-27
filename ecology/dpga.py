"""
Dual-Population Genetic Algorithm (dgpa) for Adaptive Diversity Control
"""

m = 10 # main population initial amount
n = 3 # reserve population initial amount

M = [m] # list of main population amounts at each time interval
R = [n] # same for reserve pop

def dpga():
    """
    Maintain population diversity with genetic algorithms. Dual-population
    genetic algorithm does that using an extra population for more diversity.
    """
    
