"""
Semi-global (semi global semiglobal) alignment smgb.
"""

def smgb(v, w, sigma):
    """
    Returns the semiglobal alignment of v and w and the associated score.
    """
    # Initialize the matrices.
    S = [[0 for j in xrange(len(w)+1)] for i in xrange(len(v)+1)]
    backtrack = [[0 for j in xrange(len(w)+1)] for i in xrange(len(v)+1)]
