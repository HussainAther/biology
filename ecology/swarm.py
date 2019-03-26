
"""
In swarm intelligence and computational intelligence and metaheuristics, the Bees Algorithm
relates the optimal colony formation for the foraging behavior of honey bees. As honey bees 
collecte nector from vast areas around their hive, bees communicate via a waggle dance
to inform other bees to the direction, distance, and quality rating of food sources.

We use information processing to locate and explore good sites in a problem search space.
We randomly sample the problem space to locate good sites via a local search in which
a small number of good sites are explored more frequently than others. We continue to 
exploit good sites.
"""

def evalsum(s, x):
    """
    Evaluate sum s and element x in a vector.
    """
    return s + (x**2)

def objfunc(v):
    """
    Some objective function to act upon a vector. We use this to determine
    how well our algorithm is functioning.
    """
    return all(evalsum(s, x) for s, x in v)
