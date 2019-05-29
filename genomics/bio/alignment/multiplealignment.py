"""
Multiple alignment.
"""

def multiple_alignment(word_list):
    """
    Return the multiple alignment of a given list of words.
    """
    from itertools import product
    from operator import add, mul
    from scipy.misc import comb

