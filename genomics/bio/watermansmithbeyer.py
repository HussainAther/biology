import numpy as np

"""
The dynamic programming approach by Michael S. Waterman, Temple F. Smith and William A. Beyer (1976)
computes optimal global alignments for two sequences and allows an arbitrary scoring of consecutive
gaps (insertions/deletions). It can incorporate the assumption that a single large insertion/deletion
event is biologically more likely to happen compared to many small insertions/deletions by using e.g. a
logarithmic or affine gap scoring function.

To this end, all possible gap lengths are explicitly considered within the recursion, which increases the time complexity from
O(n^2) to O(n^3) when compared to the approach by Needleman and Wunsch (1970) that applies a linear gap scoring. Under the assumption that both input sequences
a and b stem from the same origin, a global alignment tries to identify matching parts and the changes needed to transfer one sequence into the other.
The changes are scored and an optimal set of changes is identified, which defines an alignment. The dynamic programming approach
tabularizes optimal subsolutions in matrix D. An entry D_ij represents the best score for the alignment of the prefixes a_1i with b_ij. We use adding recursion.


"""
def loopstep(M, b, bb='', old_i=0):
    """
    Flip the matrix M to get the last occurence of the maximum value of M.
    Start at the highest scoring matrix cell and go until a cell with score zero
    is encountered, which should be the highest scoring local alignment.
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

def wsb(a, b, match_score=3, gap_cost=2):
    """
    Use itertools instead of two different arrays. Optimize similarity measure locally first, then
    perform the loopstep (traceback) to get the highest scoring local alignment.
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
