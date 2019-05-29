"""
Local alignment with affine gap penalty (LAFF).
"""

def local_alignment_affine_gap_penalty(v, w, scoring_matrix, sigma, epsilon):
    """
    Return the score and local alignment substrings for strings v, w with the
    given scoring matrix, gap opening penalty sigma, and gap extension penalty epsilon.
    """
