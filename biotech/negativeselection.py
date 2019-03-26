from random import random

"""
The Negative Selection algorithm is based off the self-nonself discrimintaino behavior
in the mammalian acquired immune system.
"""

def randvec(minmax):
    """
    Return a random vector from minmax.
    """
    result = []
    for i in minmax:
        result.append(minmax[i][0] + ((minmax[i][1] - minmax[i][0]) * random())) 
    return result

def eucdist(a, b):
    """
    Return Euclidean distance between two vectors a and b.
    """

