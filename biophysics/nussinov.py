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

def optimal(i, j, seq):
    """
    Return optimal pairing score between two indices i and j.
    """
    if i >= j - 4: # no pairs when i and j are fewer than four bases apart
        return 0

    
