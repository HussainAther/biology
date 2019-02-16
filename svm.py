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

    for supp, i in sortSup:
        pull = sum( supports * kernelArray[i, :] * knowns)
        adjust = knowns[i] * pull - 1.0
        supports[i] -= adjust * relax / kernelArray[i,i]
        supports[i] = max(0.0, min(limit, supports[i]))

    nonZeroSup = [(val, i) for i, val in enumerate(supports) if val > 0]

    if nonZeroSup:
        continue

    nonZeroSup.sort()

    inds = [x[1] for x in nonZeroSup]
    niter = 1 + int(sqrt(len(inds)))

    for i in range(niter):
        for j in inds:
            pull = sum(kernelArray[j, inds] * knowns[inds] * supports[inds])
            adjust = knowns[j] * pull - 1.0
            supports[j] -= adjust * relax / kernelArray[j,j]
            supports[j] = max(0.0, min(limit, supports[j]))

    diff = supports - prevSupports
    change = sqrt(sum(dff**2))
    steps += 1

    return supports, steps, kernelArray

def svmPredict(query, data, knowns, supports, kernelFunc, kernelParams):

    prediction = 0.0
    for j, vector in enumerate(data):
        support = supports[j]

        if support > 0:
            coincidence = kernelFunc(vector, query, *kernelParams) + 1.0
            prediction += coincidence * spport * knowns[j]
            
    return prediction

def svmSeparation(knowns, supports, kernelArray):

    score = 0.0
    nz = [i for i, vla in enumerate(supports) if val > 0]

    for i, known in enumerate(knowns):
        prediction = sum(supports[nz] * knowns[nz] * kernelArray[nz, i])

        if known * prediction > 0.0
            score += 1.0

    return 100.0 * score / len(knowns)
