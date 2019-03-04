import numpy as np

"""
Motoo Kimura's model of population genetics uses more than just the alpha parameter (which indicates the
general probability of transversion). It uses beta as the probability of inversions. This means the
substitution probability is beta + beta + alpha (two possible transversions and one possible substitution).
"""

def geneticDistance(p, q):
    """
    Estimate the genetic distance based on the Jukes-Cantor model updated with the Kimura model.
    p and q are respectively the fraction of transitions and transversions.
    """
    return -(1/2)*np.log(1 - 2*p - 1) - (1/4)*np.log(1 - 2*q)
