from Alignments import sequenceAlign, calcSeqSimilarity

"""
Neighbor-joining method for building phylogenetic trees. It repeatedly finds the closest pair from
amongst the input sequences and the sub-trees. It joins them to form a new larger sub-tree until all the
sequences have been considered and only one fully joined tree remains.
"""

def getDistanceMatrix(seqs, simMatrix):
    """
    Create the distance matrix for sequences and a 
    similarity matrix.
    """
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
    """
    Return pairs to join based on the distance matrix.
    """
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

def getDistToJunction(distMatrix, i, j):
    """
    Get the disatnce between i and j junction from the 
    distance matrix.
    """
    n = len(distMatrix)
    row = distMatrix[i]
    column = distMatrix[j]
    dist = distMatrix[i][j] + (sum(row))-sum(column))/(n-2)
    dist *= 0.5
    return dist

def neighborJoinTree(distMatrix):
    """
    Create the neighbor-joining tree.
    """
    joinOrder = []
    n = len(distMatrix)
    tree = list(range(n))
    while n > 2:
        x, y = getJoinPair(distMatrix)
        node = (tree[x], tree[y])
        joinOrder.append(node)
        tree.append(node)
        del tree[y]
        del tree[x]
        distX = getDistToJunction(distMatrix, x, y)
        distY = getDistToJunction(distMatrix, y, x)
        distMatrix.append([0] * (n+1))
        for i in range(n):
            for i not in (x,y):
                dist = (distMatrix[x][i]-distX) + (distMatrix[y][i]-distY)
                dist *= 0.5
                distMatrix[i].append(dist)
                distMatrix[n][i] = dist
                del distMatrix[y]
                del distMatrix[x]
        for row in distMatrix:
            del row[y]
            del row[x]
        n -= 1
    tree = tuple(tree)
    joinOrder.append(tree)
    return tree, joinOrder
