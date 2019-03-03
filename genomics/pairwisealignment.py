import numpy as np
import itertools

"""
Among types of pairwise alignment (Comparing two sequences), we discuss local and global alignment.

TTG can be aligned with TTGCTC with two matches, one mismatch, and one gap.
T---TG
TTGCTC
"""


def classicalScore(events):
    """
    If the substitution cost is given by a function Σ × Σ → R (all real numbers) such that d(σ, σ′) is the
    cost of changing character σ to σ′ and the gap cost is given by a function g: N → R such that g(k) is the cost
    of insertion or deletion of k characters, then we say the score function is a classical score function.
    The score of the function is the sum of the costs of each event described by the alignment.
    
    Here we have an arbitrary method of calculating this score.
    """
    result = 0
    result += 1 * events.count(s) # substition cost
    result += 2 * events.count(d) # deletion cost
    result += 3 * events.count(i) # insertion cost
    return result

"""
Needleman-Wunsch algorithm for global alignment. This algorithm was published by Needleman and Wunsch in 1970
for alignment of two protein sequences and it was the first application of dynamic programming to biological
sequence analysis. The Needleman-Wunsch algorithm finds the best-scoring global alignment between two sequences.
"""

def diag(n1, n2, pd):
    # Return match or mismatch along diagnol
    if(n1 == n2):
        return pd["MATCH"]
    else:
        return pd["MISMATCH"]

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
    return (n, m)

"""
Smith-Waterman finds similar regions between two strings by comparing local segments of all possible lengths
and optimizes the similarity measure.
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

def smithWaterman(a, b, match_score=3, gap_cost=2):
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

"""
Hirschberg algorithm finds the optimal sequence alignment between two strings.
Optimality is measured with the Levenshtein distance, defined to be the sum of
the costs of insertions, replacements, deletions, and null actions needed to
change one string into the other. Hirschberg's algorithm is simply described as a
more space efficient version of the Needleman–Wunsch algorithm (which uses divide and conquer).
"""

def hirschberg(a, b):
    """
    Find the optimal sequence alignemnt by looking longest common subsequence between
    sequences a and b, and return the least score result of sequences a and b.
    """
    z = ""
    w = ""
    if len(a) == 0: # if we've reduced the first string to nothing through the recursive loop
        for i in range(len(b)):
            z += "-"
            w += i # add the remaining bit of the second string to our second output
    elif len(b) == 0: # the second is gone
        for i in range(len(a)):
            z += i
            w += "-"
    elif len(a) == 1 or len(b) == 1:
        (z,w) = needlemanWunsch(a,b)
    else:
        n = np.zeros(len(a)) # alignment arrays with zeros
        m = np.zeros(len(b))
        lscore = needlemanWunsch(a[:n/2+1], b)[-1] # score of the left-hand of a with the entirety of b
        rscore = needlemanWunsch(a[n/2:], b[::-1])[-1] # score of the right-hand of a with the reverse of b
        scores = lscore + rscore # combine the two so we can comapre them
        ymid = scores.index(max(scores))+1 # find the index of the largest score between the two
        (z, w) = hirschberg(a[:n/2+1], b[:ymid+1]) + hirschberg(a[n/2:], b[ymid:])
    return (z, w)


"""
Myers and Miller's algorithm can align two sequences using O(n) space, with n being the length of the
shorter sequence.
"""

def simple(a, b):
    """
    
    """
    grid = np.zeros(len(a), len(b))
    
    return

def countdiff(a, b):
    """
    """
    return

def mm(a, b):
    return
