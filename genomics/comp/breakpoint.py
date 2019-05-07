"""
We may use breakpoints when dealing with permutations to figure out
what the identity permutation (ancestor) from which they originated.
We assume the ancestor doesn't have any breakpoints and note the 
breakpoints that appear as reversals occur in the sequence. 
"""

def brs(a):
    """
    For an array a of sequences, remove as many breakpoints as possible
    in every step.
    """
