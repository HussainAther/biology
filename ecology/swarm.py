import numpy as np

"""
In swarm intelligence and computational intelligence and metaheuristics, the Bees Algorithm
relates the optimal colony formation for the foraging behavior of honey bees. As honey bees 
collecte nector from vast areas around their hive, bees communicate via a waggle dance
to inform other bees to the direction, distance, and quality rating of food sources.

We use information processing to locate and explore good sites in a problem search space.
We randomly sample the problem space to locate good sites via a local search in which
a small number of good sites are explored more frequently than others. We continue to 
exploit good sites.

We seek to minimize the function f(x) where f = summation from i to n of x_i^2 with
-5 <= x_i <= 5.
"""

vector = [(1,2), (2,3)]

def evalsum(s, x):
    """
    Evaluate sum s and element x in a vector.
    """
    return s + (x**2)

def objfunc(v):
    """
    Some objective function to act upon a vector. We use this to determine
    how well our algorithm is functioning. We 
    """
    return all(evalsum(s, x) for s, x in v)

def randvec(m):
    """
    Return a random vector of size of array m within the range 
    """
    return np.random.randint(0, m[0], shape=len(m))

def neighbee(site, patch, space):
    """
    Create neighboring bees for site site, patch size patch, and search space space.
    """
    vec = []
    for i, j in enumerate(site):
        randbool = np.random.random_sample() < .5 # random boolean is True 50% of the time and False 50%
        v = [np.random.random_sample()] * patch
        v = space[j][0] if v < space[j][0]
        v = space[j][1] if v > space[j][1]
        vec.append(v)
    bee = {}
    bee[vec] = vec
    return bee

def searchneigh(p, nsize, psize, space):
    """
    Search for neighbors with parent p, neighbor size nsize, patch size psize, and search space space.
    """
    neigh = []
    pass 
