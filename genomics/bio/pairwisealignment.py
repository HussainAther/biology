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

def diag(x, y, pd):
    # Return match or mismatch along diagnol
    if(x == y):
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
    simple method of reducing the lengths of a and b
    """
    grid = np.zeros(len(a), len(b))
    for i in range(1, len(b) + 1):
        grid[i][0] = grid[i-1][0] + 1 # move it and increase the number
        for j in range(1, len(a) + 1):
            grid[i][j] = min(grid[i-1][j] + 1, grid[i][j-1] + 1, grid[i-1][j-1])
            if a[j-1] == b[i-1]:
                grid[i][j] = min(grid[i][j], grid[i-1][j-1])
    result = []
    i, j = len(a), len(b)
    while i != 0 and j != 0: # loop through this till we cover the lengths of both a and b.
        if i == 0:
            result.append(("+", b[j-1]))
            j -= 1
        elif j == 0:
            result.append(("-", a[i-1]))
            i -= 1
        else:
            diag = grid[j-1][i-1]
            up = grid[j-1][i] + 1
            left = grid[j][i-1] + 1
            curr = grid[j][i]
            if curr = diag and a[i-1] == b[j-1]:
                result.append((" ", b[i-1]))
                j -= 1
                i -= 1
            elif: curr = left:
                result.append(("-", b[i-1]))
                i -= 1
            else:
                result.append(("+", b[j-1]))
                j -= 1
    return result[::-1]

def countdiff(a, b):
    """
    reduce the lengths through the counting method. compare two sequences row and rw against one another
    and return the last row to show the differences between the sequences
    """
    row = range(0, len(a)+1)
    for i in range(1, len(b)+1):
        rw = [] * (len(a)+1)
        rw[0] = i
        for j in range(1, len(a) + 1):
            rw[j] = min(rw[j-1] + 1, row[j] + 1)
            if a[j-1] == b[i-1]:
                rw[j] = min(rw[i], row(j-1]))
        row, rw = rw, row # switch them till we get the last row
    return row

def mm(a, b):
    """
    Miller-Myers algorithm is a recursive algorithm based on the Hirschberg method that
    exploits the lengths of the sequences to perform at a greater speed..
    """
    if len(a) < 3 or len(b) < 3:
        return simple(a, b)
    if b % 2 == 0: # find the midpoint mid
        mid = len(b)/2
    else:
        mid = len(b)/2 + 1
    row1 = countdiff(a, b[:mid]) # first row
    st1 = a[::-1] # reverse the first sequence
    st2 = b[(len(b) - mid+1):][::-1] # reverse the latter half of b
    row2 = countdiff(st1, st2)[::-1] # count the differences between them and reverse the output row
    index = 0 #
    for i in range(1, len(a)+1):
        if row1[index] + row2[index] > row1[i] + row2[i]: # find an index at which the diffrences at that index are greatest
            index = i
    return mm(a[:index], b[:mid]) + mm(a[index:], b[mid:]) # loop again

"""
The Method of Four Russians partitions the matrix into small square blocks and uses a lookup
table to perform the algorithm quickly within each block. The lookup table encodes an index from the upper
left quadrant of the block to the boundary cells on the lower right of the block. It improves the overall
speed from (n/t)^2 blocks instead of n^2 matrix cells in which n is the side length of the matrix.
"""

lookup = {}

def fourRussians(a, b):
    return

"""
The Gotoh algorithm computes the optimal global alignment of two sequences when using an affine gap scoring.
Here, the scoring of a long consecutive gap (insertion/deletion) is favored over a collection of small gaps
with the same combined length. This incorporates the assumption that a single large insertion/deletion event
is biologically more likely to happen compared to many small insertions/deletions. While sophisticated gap
scoring models can be applied in the generic algorithm by Waterman-Smith-Beyer (1976), affine gap scoring used
in Gotoh's algorithm enables a reasonable gap model with reduced runtime.
"""

def affine(i, j):
    """
    Affine scoring for two indices i and j.
    """
    cons = 5 # some constant
    return i + j * cons

def gotoh(a, b):
    """
    Use itertools instead of two different arrays. Optimize similarity measure locally first, then
    perform the loopstep (traceback) to get the highest scoring local alignment. Use affine scoring to
    measure penalties
    """
    a, b = a.upper(), b.upper() # convert to uppercase
    M = np.zeros((len(a)+1, len(b)+1))
    prev = "" # previous indel
    for i, j in itertools.product(range(1, M.shape[0]), range(1, M.shape[1])):
        c = M[i - 1, j - 1] + (match_score if a[i - 1] == b[j - 1] else - match_score) # score if there's a match
        d = M[i - 1, j] - gap_cost # deletion
        v = M[i, j - 1] - gap_cost # insertion
        M[i, j] = max(c, d, v, 0) # add the maximum of the possible scores
        if M[i, j] == prev: # if the indel is the same as the one before
            blockgap_cost = affine(i, j) # affine scoring
            M[i, j] == blockgap_cost # use the block gap cost instead
        prev = max(c, d, v, 0) # get the previous indel
    bb, pos = loopstep(M, b) # loop through the second string locally until the final matrix index is 0
    return (pos, pos + len(bb))


"""
Waterman-Smith-Beyer algorithm.

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
    blockgap_cost = 1 # if there's a block of gaps, then that's worth a leser cost than several individual gaps
    a, b = a.upper(), b.upper() # convert to uppercase
    M = np.zeros((len(a)+1, len(b)+1))
    prev = "" # previous indel
    for i, j in itertools.product(range(1, M.shape[0]), range(1, M.shape[1])):
        c = M[i - 1, j - 1] + (match_score if a[i - 1] == b[j - 1] else - match_score) # score if there's a match
        d = M[i - 1, j] - gap_cost # deletion
        v = M[i, j - 1] - gap_cost # insertion
        M[i, j] = max(c, d, v, 0) # add the maximum of the possible scores
        if M[i, j] == prev: # if the indel is the same as the one before
            M[i, j] == blockgap_cost # use the block gap cost instead
        prev = max(c, d, v, 0) # get the previous indel
    bb, pos = loopstep(M, b) # loop through the second string locally until the final matrix index is 0
    return (pos, pos + len(bb))

"""
Feng-Doolittle algorithm.

The progressive alignment approach by Da-Fei Feng and Russell F. Doolittle (1987) computes a multi-sequence-alignment (MSA)
of a set of sequences based on pairwise alignments. This approximative approach identifies good MSA solutions in reasonable time.

"""

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

def fd(a, b):
    """
    Compute pairwise alignments of all sequences, use similarities instead of distances,
    and build a tree with the values. Then use the tree to generate alignments.
    """
    a, b = a.upper(), b.upper() # convert to uppercase
    M = np.zeros((len(a)+1, len(b)+1))
    for i, j in itertools.product(range(1, M.shape[0]), range(1, M.shape[1])):
        c = M[i - 1, j - 1] + (match_score if a[i - 1] == b[j - 1] else - match_score) # score if there's a match
        d = M[i - 1, j] + penaltydict["GAP"]# deletion
        v = M[i, j - 1] + penaltydict["GAP"] # insertion
        M[i, j] = max(c, d, v, 0) # add the maximum of the possible scores
    bb, pos = loopstep(M, b) # loop through the second string locally until the final matrix index is 0
    return (pos, pos + len(bb))
