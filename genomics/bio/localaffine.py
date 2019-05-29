"""
Local alignment with affine gap penalty (LAFF).
"""

def local_alignment_affine_gap_penalty(v, w, scoring_matrix, sigma, epsilon):
    """
    Return the score and local alignment substrings for strings v, w with the
    given scoring matrix, gap opening penalty sigma, and gap extension penalty epsilon.
    """
    # Initialize the matrices
    S_lower = [[0 for j in xrange(len(w)+1)] for i in xrange(len(v)+1)]
    S_middle = [[0 for j in xrange(len(w)+1)] for i in xrange(len(v)+1)]
    S_upper = [[0 for j in xrange(len(w)+1)] for i in xrange(len(v)+1)]
    backtrack = [[0 for j in xrange(len(w)+1)] for i in xrange(len(v)+1)]

    # Initialize the maximum score below the lowest possible score.
    max_score = -1
    max_i, max_j = 0, 0

    # Fill in the Score and Backtrack matrices.
    for i in xrange(1, len(v)+1):
        for j in xrange(1, len(w)+1):
            S_lower[i][j] = max([S_lower[i-1][j] - epsilon, S_middle[i-1][j] - sigma])
            S_upper[i][j] = max([S_upper[i][j-1] - epsilon, S_middle[i][j-1] - sigma])
            middle_scores = [S_lower[i][j], S_middle[i-1][j-1] + scoring_matrix[v[i-1], w[j-1]], S_upper[i][j], 0]
            S_middle[i][j] = max(middle_scores)
            backtrack[i][j] = middle_scores.index(S_middle[i][j])
