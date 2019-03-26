import numpy as np

"""
We use the Particle Swarm intelligence (particle swarm, PSO) algorithm to simulate behaviors
of swarms and optimize these behaviors iteratively. We search for a best solution by moving
particles with a velocity and calculating the solution every iteration until the particle swarm
converges to the best solution.
"""
 
x = [] # x-value positions
v = [] # velocity values

def fit(i):
    """
    Return fitness value for some particle index i .
    """
    return x[i] ** 2 + 1

def pso(a, dim):
    """
    For an array of particles a with dimension dim, we search for a solution that maximizes  
    the fitness value for each particle.
    """
    x = np.random.rand(len(a))*100 # x positions of each particle
    v = np.random.rand(len(a))*100 # v velocity of each particle
    p_best = np.zeros(len(a)) # best fitness values of each particle
    k = 1 # iteration
    for i,j in enumerate(a):
         if fit(i) > p_best[i]:
             p_best[i] = fit(i) # calculate the current fitness values
    g_best = p_best.index(max(p_best)) # the particle that is the most fit

