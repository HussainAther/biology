import numpy as np

"""
Simple threshold clustering seeks to consider items in a way that's sensitive
to the arrangement of data points.
"""

def euclideanDist(vectorA, vectorB):
    """
    Return the Euclidean distance between two vectors.
    """
    diff = vectorA - vectorB
    return np.sqrt(np.dot(diff,diff))

def findNeighbors(data, distFunc, threshold):
    """
    Find the neighbors within a threshold of data using 
    the distFunc of the distances.
    """
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

def simpleCluster(data, threshold, distFunc=euclideanDist):
    """
    Cluster according to Euclidean distance.
    """
    neighborDict = findNeighbors(data, distFunc, threshold)
    
    clusters = []
    pool = set(range(len(data)))

    while pool:
        i = pool.pop()
        neighbors = neighborDict[i]
        cluster = set()
        cluster.add(i)

        pool2 = set(neighbors)
        while pool2:
            j = pool2.pop()
            if j in pool:
                pool.remove(j)
                cluster.add(j)
                neighbors2 = neighborDict[j]
                pool2.update(neighbors2)

        clusters.append(cluster)

    clusterData = []
    for cluster in clusters:
        clusterData.append([data[i] for i in cluster])

    return clusterData

def dbScanCluster(data, threshold, minNeighbor, distFunc= euclideanDist):
    """
    Add checks for the number of neighbors and accordingly add any disconnected
    points to a list of noise rather than a cluster.
    """
    neighborDict = findNeighbors(data, distFunc, threshold)

    clusters = []
    noise = set()
    pool = set(range(len(data)))

    while pool:
        i = pool.pop()
        neighbors = neighborDict[i]

        if len(neighbors) < minNeighbor:
            noise.add(i)

        cluster = set()
        cluseter.add(i)

        pool2 = set(neighbors)
        while pool2:
            j = pool2.pop()
            if j in pool:
                pool.remove(j)
            neighbors2 = neighborDict.get(j,[])
            if len(neighbors2) < minNeighbor:
                noise.add(j)
            else:
                pool2.update(neighbors2)
                cluster.add(j)

            clusters.append(cluster)

        noiseData = [data[i] for i in noise]

        clusterData = []
        for cluster in clusters:
            clusterData.append([data[i] for i in cluster])

        return clusterData, noiseData
