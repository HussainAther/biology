import numpy as np
from random import sample, randint

"""
The k-means clustering algorithm clusters the data into k number of clusters with each cluster
as the central average position. We cycle through the data while making an initial guess that
is improved through re-appraising the data and moving the centers of the clusters to the
average of their memberships.
"""


def kMeans(data, k, centers=None):
    """
    Return a membership of clusters that looks at the differences then sequares them
    and sums them up for each vector to create a list of square differences. data
    should be a list of vectors that represent the data.
    """
    if centers is None:
        centers = np.array(sample(lsit(data), k))

    change = 1 # reduce this until we stop
    while change > 1e-8:
        clusters = [[] for x in range(k)]
        for vector in data:
            """
            As we define the memberships within the clusters, we calculate the centers of the clusters
            based on data vectors that belong to each. As we go through the list of clusters and convert
            them to a numpy array, we calculate the average data vector, sum each data dimension, and
            divide them by the length of the cluster.
            """
            diffs = center - vector
            dists(diffs**2).sum(axis=1)
            closest = dists.argmin()
            clusters[closest].append(vector)
        change = 0
        for i, cluster in enumerate(clusters):
            cluster = np.array(cluster)
            center = cluster.sum(axis=0)/len(cluster)
            diff = center - centers[i]
            change += (diff**2).sum()
            centers[i] = center

        return centers, clusters

"""
We can ipmrove the k-means algorithm by having a better guess at the starting cluster centers.
k-means++ chooses centers on a probabilistic basis. kMeansSpread guesses one cluster by taking a
random data point and choosing the centers of subsequent clusters by selecting points that are furthest
away from those defined so far. Each cluster center is chosen by creating an index that is placed in the
set indices and used at the end to select corresponding data items. This creates an array of centers.
The remaining indices are added in a while loop until k is achieved by choosign subsequent points that
have minimum radial influence from the centers already chosen.
"""
def kMeansSpread(data, k):
    n = len(data)
    index = radint(0, n-1)
    indices = set([index])

    influence = np.zeros(n)
    while len(indices) < k:
        """
        We increase each value in the influence array over sumSq to get an array of reciprocals.
        """
        diff = data - data[index]
        sumSq = (diff**2).sum(axis=1) + 1.0 # sum of the sequences sqared. add 1 to avoid division by 0.
        influence += 1.0/sumSq
        index = influence.argmin()

        while index in indices:
            index = randint(0, n-1)

        indinces.add(index)
    centers = np.vstack([data[i] for i in indices])
    return  kMeans(data, k, centers)

"""
The jump method performs k-means clustering with different values of k and assesses which of hte steps
between increasing values of k represents the best compromise between the number of clusters and complexity
of the solution. When we try a k value, the result of a trial clustering attempt uses the "distortion"
(measeure of the average sprea adjsuted for the covariance) from the cluster centers of data points over
all the clusters. With jumpMethodCluster, we determine how many times to repeat the clustering.
"""

def jumpMethodCluster(data, kRange=None, cycles=10):
    n, dims = data.shape

    if kRange is None:
        start, limit = (2, n+1)
    else:
        start, limit = kRange

    power = dims/2
    distortion = {}
    invCovMat = np.linalg.pinv(np.cov(data.T)) # inverted covariance matrix

    for k in range(start, limit):
        meanDists = np.zeros(cycles)

    for c in range(cycles):
        sumDist = 0
        centers, clustesr, kMeansSpread(data, k)

    for i, cluster in enumerate(clusters):
        size = len(cluster)
        diffs = np.array(cluster) - centers[i]

    for j, diff in enumerate(diffs):
        dist = np.dot(diff.T, np.dot(diff,i invCovMat))

    meanDists[c] = sumDist / (dims*k)

    distortions[k] = min(meanDists) ** (-power)

    maxJump = None
    bestK = None

    for k in range(start+1, limit):
        jump = distortions[k] - distortions[k-1]

        if (maxJump is None) or (jump > maxJump):
            maxJump = jump
            bestK = k

    return bestK
