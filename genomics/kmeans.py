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
            diffs = center - vector
            dists(diffs**2).sum(axis=1)
            closest = dists.argmin()
            clusters[closest].append(vector)
