import numpy as np

"""
We use the Particle Swarm intelligence (particle swarm, PSO) algorithm to simulate behaviors
of swarms and optimize these behaviors iteratively. We search for a best solution by moving
particles with a velocity and calculating the solution every iteration until the particle swarm
converges to the best solution.
"""
 
p = [] # positions
v = [] # velocity values

def fit(i, dim):
    """
    Return fitness value for some particle index i and dimensions dim.
    """
    result = 0
    for j in range(dim):
        result += p[j][i] ** 2 + 1
    return result

def pso(a, dim):
    """
    For a number of particles with dimension dim, we search for a solution that maximizes  
    the fitness value for each particle.
    """
    for i in range(dim):
        p.append(np.random.rand(a)*100) # positions of each particle
    for i in range(dim):
        v.append(np.random.rand(x)*100) # each velocity component for each particle
    p_best = np.zeros(a) # best fitness values of each particle
    for i in range(a):
        if fit(i, dim) > p_best[i]:
            p_best[i] = fit(i, dim) # calculate the current fitness values
    g_best = p_best.index(max(p_best)) # the particle that is the most fit
    for i in range(a):
        for j in range(dim):
            v[i] += np.random.rand(p[i][j]-p[1][j]) + np.random.rand(p[g_best][j] - p[1][j]) # update velocity
            p[i][j] += v[i] # update position
    return p, v

