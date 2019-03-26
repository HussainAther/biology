"""
We use the Particle Swarm intelligence (particle swarm, PSO) algorithm to simulate behaviors
of swarms and optimize these behaviors iteratively. We search for a best solution by moving
particles with a velocity and calculating the solution every iteration until the particle swarm
converges to the best solution.
"""

def pso(a, dim):
    """
    For an array of particles a with dimension dim, we search for a solution that maximizes  
    the fitness value for each particle.
    """
