from itertools import product
from operator import add, mul
from scipy.misc import comb

"""
Multiple alignment.
"""

def multiple_alignment(word_list):
    """
    Return the multiple alignment of a given list of words.
    """
    word_list = ["$" + word for word in word_list]

    # Initialize scoring and backtrack dictionaries, along with the indices and base score.
    S, backtrack = {}, {}
    perm_list = list(product([0, -1], repeat=len(word_list)))[1:]
    base_score = -1*comb(len(word_list), 2, exact=True)

