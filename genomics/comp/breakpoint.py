"""
We may use breakpoints when dealing with permutations to figure out
what the identity permutation (ancestor) from which they originated.
We assume the ancestor doesn't have any breakpoints and note the 
breakpoints that appear as reversals occur in the sequence. 
"""

def issorted(l):
    """
    Is the list l sorted?
    """
    return sorted(l) == list(range(min(l), max(l)+1))

def brs(a):
    """
    For an array a of sequences of numbers, remove as many breakpoints as possible
    in every step until we get the numbers in order. 
    """
    for i in a:
        while not issorted(i):
 
