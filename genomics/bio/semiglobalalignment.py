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

    # Fill in the Score and Backtrack matrices.
    for i in xrange(1, len(v)+1):
        for j in xrange(1, len(w)+1):
            scores = [S[i-1][j] - sigma, S[i][j-1] - sigma, S[i-1][j-1] + [-1, 1][v[i-1] == w[j-1]]]
            S[i][j] = max(scores)
            backtrack[i][j] = scores.index(S[i][j])

    # Get the position of the highest scoring cell in the last row or last column.
    last_row_index = max(xrange(len(w)+1), key=lambda x: S[len(v)][x])
    last_column_index = max(xrange(len(v)+1), key=lambda x: S[x][len(w)])
    if S[len(v)][last_row_index] >= S[last_column_index][len(w)]:
        i = len(v)
        j = last_row_index
    else:
        i = last_column_index
        j = len(w)
    max_score = S[i][j]

