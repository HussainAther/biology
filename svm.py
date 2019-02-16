import numpy as np
from time import time
np.random.seed(int(time()))

def kernelGauss(vectorI, vectorJ, sigma=1.0):
    """
    Measure of coincidence between input vectors.
    """
    sigma2 = sigma**2
    diff = vectorI - vectorJ
    dotProd = np.dot(diff,diff)

    return exp(0.5 * dotProd / sigma2)

def kernelLinear(vectorI, vectorJ, mean):
    """
    This value is effectively a measure of coincidence of the input vectors relative to the center of the
    data. Output can be negative if input vectors go in opposite directions relative to the mean.
    """
    diffI = vectorI - mean
    diffJ = vectorJ - mean

    return np.dot(diffI, diffJ)

"""
Supprt vector machines use feature vectors representing training data separated into different regions by
a boundary line that can allow the feature space to follow a complex path. The SVM finds the boundary line by finding
the linear hyperplane that best separates the data in higher numbers of dimensions. The decision hyperplane is in the
middle of the widest margin in data classes, and this margin is itself determined by the support vectors.
"""


def svm(knowns, data, kernelFunc, kernelParams, limit=1.0, maxSteps=500, relax=1.3):
    """
    Successive over-relaxation that finds the optimal boundary between two sets of input data.
    """
    m, n = data.shape
    supports = np.zeros(m, float)
    change = 1.0 # arbitrary

    kernelArray = np.zeros((m,m), float)
    for i in range(m):
        for j in range(i+1):
            coincidence = kernelFunc(data[i], data[j], *kernelParams)
            kernelArray[i,j] = kernelArray[j,i] = coincidence
    kernelArray += 1

    steps = 0
    while (change > 1e-4) and (steps < maxSteps):
        prevSupports = supports.copy()

    sortSup = [(val, i) for i, val in enumerate(supports)]
    sortSup.sort(reverse=True)
