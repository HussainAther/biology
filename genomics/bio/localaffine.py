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
