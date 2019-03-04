import numpy as np

"""
Simple threshold clustering seeks to consider items in a way that's sensitive
to the arrangement of data points.
"""

def euclideanDist(vectorA, vectorB):
    diff = vectorA - vectorB
    return np.sqrt(np.dot(diff,diff))

def findNeighbors(data, distFunc, threshold):
    neighborDict = {}
    n = len(data)
    for i in range(n):
        neighborDict[i] = []

    for i in range(0, n-1):
        for j in range(i+1, n):
            dist = distFunc(data[i], data[j])

            if dist < threshold:
                neighborDict[i].append(j)
                neighborDict[j].append(i)

    return neighborDict
