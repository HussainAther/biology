import numpy as np

"""
Local (local) alignment

Initialize first row and first column to be 0. The score of the best local alignment is the largest value
in the entire array. To find the actual local alignment: start at an entry with the maximum score, traceback 
as usual, stop when we reach an entry with a score of .0
"""

def local_align(x, y, score=ScoreParam(-7, 10, -5)):
    """
    Do a local alignment between x and y.
    """
    A = np.zeros((len(x) + 1, len(y) + 1)) # create a zero-filled matrix
    best = 0
    optloc = (0,0)
    # fill in A in
    for i in xrange(1, len(y)):
        for j in xrange(1, len(x)):
            # the local alignment recurrance rule:
            A[i][j] = max(
            A[i][j-1] + score.gap,
            A[i-1][j] + score.gap,
            A[i-1][j-1] + (score.match if x[i] == y[j] else score.mismatch), 0)
            # track the cell with the largest score
            if A[i][j] >= best:
                best = A[i][j]
                optloc = (i,j) 
    # return the opt score and the best location
    return best, optloc
