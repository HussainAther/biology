import numpy as np

"""
A tandem repeat Pk of a pattern P = p1 . . . pm is a pattern of length k ·m formed
by concatenating k copies of P. The tandem repeat problem is to find the best
global alignment between some interval of a text T = t1 . . . tn and a tandem repeat
Pk for some k. When computing global alignments the premium for matching is
+1 and the mismatch and gap penalties are −1.

For example given the pattern P = GGT and the text T = ACCGGTGCTGTAA the best
tandem repeat is given by k = 3 and the following alignment:

ACCGGTGCTG-TAA
   GGTGGTGGT
"""

class ScoreParam:
    """
    The parameters for an alignment scoring function.
    """
    def __init__(self, gap, match, mismatch):
        self.gap = gap
        self.match = match
        self.mismatch = mismatch


def local_align(x, y, score=ScoreParam(-7, 10, -5)):
    """
    Fill an alignment matrix A with normal recursion for the tandem repeat alignment problem.
    For some pattern x and text y, we determine the best tandem repeat when aligning
    x to y. Use a local alignment between x and y.
    """
    A = np.zeros((len(x) + 1, len(y) + 1)) # create a zero-filled matrix
    best = 0
    optloc = (0,0)
    for i in xrange(1, len(y)): # fill in our matrix using the alignment
        for j in xrange(1, len(x)): # the local alignment recurrence recursion rule            
            A[i][j] = max(
            A[i][j-1] + score.gap,
            A[i-1][j] + score.gap,
            A[i-1][j-1] + (score.match if x[i] == y[j] else score.mismatch), 0)
            # track the cell with the largest score
            if A[i][j] >= best:
                best = A[i][j]
                optloc = (i,j) 
    return best, optloc # return the opt score and the best location
