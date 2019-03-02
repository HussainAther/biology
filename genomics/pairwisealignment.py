import numpy as np
import itertools

"""
Among types of pairwise alignment (Comparing two sequences), we discuss local and global alignment.

TTG can be aligned with TTGCTC with two matches, one mismatch, and one gap.
T---TG
TTGCTC
"""


def classicalScore():
    """
    If the substitution cost is given by a function Σ × Σ → R (all real numbers) such that d(σ, σ′) is the
    cost of changing character σ to σ′ and the gap cost is given by a function g: N → R such that g(k) is the cost
    of insertion or deletion of k characters, then we say the score function is a classical score function.
    The score of the function is the sum of hte costs of each event described by the alignment.
    """




"""
Needleman-Wunsch algorithm for global alignment. This algorithm was published by Needleman and Wunsch in 1970
for alignment of two protein sequences and it was the first application of dynamic programming to biological
sequence analysis. The Needleman-Wunsch algorithm finds the best-scoring global alignment between two sequences.
"""

def diag(n1, n2, pt):
    # Return match or mismatch along diagnol
    if(n1 == n2):
        return pt["MATCH"]
    else:
        return pt["MISMATCH"]

def point(d, h, v):
    # Return pointer matrix
    p = max(d, h, v)
    if(d == p):
        return "D"
    elif(h == p):
        return "H"
    else:
         return "V"


penaltydict =
    {
    "MATCH": match,
    "MISMATCH": mismatch,
    "GAP": gap
}

def needlemanWunsch(a, b, match = 1,mismatch = -1, gap = -2):
    """
    matrix is numpy array of len for rows and columns of the sizes of the two sequences.
    """
    a, b = a.upper(), b.upper() # convert to uppercase
    n = np.zeros(len(a)) # alignment arrays with zeros
    m = np.zeros(len(b))
    for i in range(len(n)): # gap penalty for first row
        n[i][0] = penaltydict["GAP"] * i
        m[i][0] = "V"
    for i in range(len(m)): # gap penalty of first column
        n[0][i] = penaltydict["GAP"] * i
        m[0][i] = "H"
    # the rest of the matrix
    for i in range(1, len(m)):
        for j in range(1, len(n)):
            d = n[i-1][j-1] + diag(a[j-1], b[i-1], penaltydict) # match or mismatch on the diagnol
            h = n[i][j-1] + penaltydict["GAP"] # gap horizontal
            v = n[i-1][j] + penaltydict["GAP"] # gap vertical
            n[i][j] = max(d, h, v) # determine which of d, h, and v have the most points
            m[i][j] = point(d, h, v) # convert it to a code
    return np.matrix(n, m)

"""
Smith-Waterman finds similar regions between two strings by comparing local segments of all possible lengths
and optimizes the similarity measure.
"""

def smithWaterman(a, b, match_score=3, gap_cost=2):
    a, b = a.upper(), b.upper() # convert to uppercase
    H = np.zeros((len(a) + 1, len(b) + 1), np.int)

    for i, j in itertools.product(range(1, H.shape[0]), range(1, H.shape[1])):
        match = H[i - 1, j - 1] + (match_score if a[i - 1] == b[j - 1] else - match_score)
        delete = H[i - 1, j] - gap_cost
        insert = H[i, j - 1] - gap_cost
        H[i, j] = max(match, delete, insert, 0)
    b_, pos = traceback(H, b)
    return pos, pos + len(b_)
