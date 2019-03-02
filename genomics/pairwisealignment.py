import numpy as np
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

penaltydict =
    {
    "MATCH": match,
    "MISMATCH": mismatch,
    "GAP": gap
}

def needlemanWunsch(matrix, match = 1,mismatch = -1, gap = -2):
    """
    matrix is numpy array of len for rows and columns of the sizes of the two sequences.
    """
    n = np.zeros(len(matrix)) # alignment matrix with zeros
    m = np.zeros(len(matrix[0]))
    for i in range(len(n)): # gap penalty for first row
        n[i][0] = penaltydict["GAP"] * i
        m[i][0] = "V"
    for i in range(len(m)): # gap penalty of first column
        n[0][i] = penaltydict["GAP"] * i
        m[0][i] = "H"
    # the rest of the matrix
    
