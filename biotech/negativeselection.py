from random import random
import numpy as np

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
    return np.linalg.norm(a-b) 
 
def contains(vec, space):
    """
    Return True if the vector is in space. False otherwise.
    """
    for i in vector:
         for j in i:
               if i< space[j][0] or i > space[j][1]:
                   return False
    return True
