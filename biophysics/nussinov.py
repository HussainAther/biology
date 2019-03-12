import sys
import numpy as np

"""
We implement that Nussinov RNA Folding Algorithm for predicting the secondary structure of RNA.
We get a two-dimensional secondary structure of RNA by determining how pairing between nucleotides
in the primary sequence give rise to a two-dimensional architecture.
"""

def pair(t):
    """
    Check if the pair of nucleotides bond with one another.
    """
    if t in [("A", "U"), ("U", "A"), ("C", "G"), ("G", "C")]:
        return True
    return False
