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

    for index in product(*map(xrange,map(lambda s: len(s) + 1, word_list))):
        # We forced a match with the first symbols, so the zero-shell should lead to the zero index.
        if reduce(mul, index) == 0:
            # Since we forced a match with the first symbol, we want to force starting point to be the zero index.
            if sum(index) == 0:
                # All symbols match.
                S[index] = 0
            else:
                # Make it smaller than the lowest possible score.
                S[index] = 2*base_score*reduce(add, map(len, word_list))
        else:
            # Use previous scores to determine the best score for the current index.
            previous_scores = [S[tuple(map(add, index, perm))] for perm in perm_list]
            current_index_scores = []
            for perm in perm_list:
                chars = [word_list[i][index[i]-1] if perm_value == -1 else '-' for i, perm_value in enumerate(perm)]
                current_index_scores.append(base_score + sum([comb(chars.count(ch), 2, exact=True) for ch in set(chars)]))
            scores = map(add, previous_scores, current_index_scores)
            backtrack[index], S[index] = max(enumerate(scores), key=lambda p: p[1])
