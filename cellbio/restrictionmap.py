from collections import Counter

"""
Create a restriction map, locations where restriction sites are in DNA.
"""

def partial_digest(distances):
    """
    Return a set whose positive pairwise differences generate distances.
    """
    # Initialize variables.
    X = {0}
    width = max(distances)

    # Create lambda functions for multiset operations.
    new_dist = lambda y, S: Counter(abs(y-s) for s in S)
    containment = lambda a, b: all(a[x] <= b[x] for x in a)
