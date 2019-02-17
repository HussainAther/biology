from Alignments import sequenceAlign, calcSeqSimilarity
from math import exp

"""
Neighbor-joining method for building phylogenetic trees. It repeatedly finds the closest pair from
amongst the input sequences and the sub-trees. It joins them to form a new larger sub-tree until all the
sequences have been considered and only one fully joined tree remains.
"""

def getDistanceMatrix(seqs, simMatrix):
    n = len(seqs)
    matrix = [[0.0] *n for x in range(n)]
    maxScores = [calcSeqSimilarity(x,x,simMatrix) for x in seqs]
