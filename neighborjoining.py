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
    for i in range(n-1):
        seqA = seqs[i]
        for j in range(i+1, n):
            seqB = seqs[j]
            score, alignA, alignB = sequenceAlign(seqA, seqB, simMatrix)
            maxScore = max(maxScores[i],maxScores[j])
            dist = maxScore - score
            matrix[i][j] = dist
            matrix[j][i] = dist
    return matrix

def getJoinPair(distMatrix):
    n = len(distMatrix)
    minQ = None
    joinPair = None
    for i in range(n-1):
        sumRow = sum(distMatrix[i])
        for j in range(i+1, n):
            sumCol = sum(distMatrix[j])
            dist = distMatrix[i][j]
            q = (n-2)*dist - sumRow - sumCol
            if (minQ is None) or (q < minQ):
                minQ = q
                joinPair = [i, j]
    return joinPair


