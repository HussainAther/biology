import numpy as np
from random import sample

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
"""
