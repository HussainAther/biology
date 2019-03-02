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

def loopstep(M, b, bb='', old_i=0):
    """
    Flip the matrix M to get the last occurence of the maximum value of M.
    """
    Mf = np.flip(np.flip(M, 0), 1) # flip it
    ii, jj = np.unravel_index(Mf.argmax(), Mf.shape) # given a linear index, find the corresponding Mf.shape-dimensinoal index
    i, j = np.subtract(M.shape, (ii+1, jj+1)) # subtract the final Mf.shape-dimensional indices from the shape of our matrix
    if M[i, j] = 0: # the final index is 0
        return bb, j # return the string and the second index
    if old_i - i > 1: # check if we ned to add a gap. if we have found ourselves one or more base ahead
        bb = b[j-1] + "-" + bb # add the gap
    else:
        bb = b[j-1] + bb
    return loopstep(M[0:1, 0:j], b, bb, i) # loop again until we find the 0 index

def smithWaterman(a, b, match_score=3, gap_cost=2):
    """
    Use itertools instead of two different arrays.
    """
    a, b = a.upper(), b.upper() # convert to uppercase
    M = np.zeros((len(a)+1, len(b)+1))
    for i, j in itertools.product(range(1, M.shape[0]), range(1, M.shape[1])):
        c = M[i - 1, j - 1] + (match_score if a[i - 1] == b[j - 1] else - match_score) # score if there's a match
        d = M[i - 1, j] - gap_cost # deletion
        v = M[i, j - 1] - gap_cost # insertion
        M[i, j] = max(c, d, v, 0) # add the maximum of the possible scores
    bb, pos = loopstep(M, b) # loop through the second string locally until the final matrix index is 0
    return (pos, pos + len(bb))
