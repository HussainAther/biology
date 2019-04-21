import numpy as np

from time import time
from matplotlib import pyplot as plt

np.random.seed(int(time()))

"""
Support vector machine (SVM svm)
"""

def kernelGauss(vectorI, vectorJ, sigma=1.0):
    """
    Measure of coincidence between input vectors.
    """
    sigma2 = sigma**2
    diff = vectorI - vectorJ
    dotProd = np.dot(diff,diff)

    return np.exp(0.5 * dotProd / sigma2)

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

def svmTrain(knowns, data, kernelFunc, kernelParams, limit=1.0, maxSteps=500, relax=1.3):
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
        pull = np.sum( supports * kernelArray[i, :] * knowns)
        adjust = knowns[i] * pull - 1.0
        supports[i] -= adjust * relax / kernelArray[i,i]
        supports[i] = max(0.0, min(limit, supports[i]))

    nonZeroSup = [(val, i) for i, val in enumerate(supports) if val > 0]

    nonZeroSup.sort()

    inds = [x[1] for x in nonZeroSup]
    niter = 1 + int(np.sqrt(len(inds)))

    for i in range(niter):
        for j in inds:
            pull = np.sum(kernelArray[j, inds] * knowns[inds] * supports[inds])
            adjust = knowns[j] * pull - 1.0
            supports[j] -= adjust * relax / kernelArray[j,j]
            supports[j] = max(0.0, min(limit, supports[j]))

    diff = supports - prevSupports
    change = np.sqrt(np.sum(dff**2))
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
        prediction = np.sum(supports[nz] * knowns[nz] * kernelArray[nz, i])

        if known * prediction > 0.0:
            score += 1.0

    return 100.0 * score / len(knowns)

numPoints = 20
catData = []

for x in range(1,6):
    for y in range(1,6):
        xNorm = x/6.0
        yNorm = y/6.0
 
        if x == 3 and y == 3: # one category or the other
            category = -1.0

        elif x % 2 == y % 2:
            category = 1.0

        else:
            category = -1.0

        xvals = np.random.normal(xNorm, 0.2, numPoints)
        yvals = np.random.normal(yNorm, 0.2, numPoints)

        for i in range(numPoints):
            catData.append((xvals[i], yvals[i], category))

catData = np.array(catData)
np.random.shuffle(catData)

knowns = catData[:, -1] 
data = catData[:, :-1]

params = (0.1,)
supports, steps, kernelArray = svmTrain(knowns, data, kernelGauss, params)

score = svmSeparation(knowns, supports, kernelArray)
print("Known data: %5.2%% correct" % (score))

ds1x = []
ds1y = []
ds2x = []
ds2y = []

x = 0.0
while x < 1.0:
    y = 0.0
    while y < 1.0:
        query = np.array((x,y))
        prediction = svmPredict(query, data, knowns, supports, kernelGauss, params)

        if prediction > 0:
            ds1x.append(x)
            ds1y.append(y)
        else:
            ds2x.append(x)
            ds2y.append(y)
        y += 0.02
    x += 0.02
plt.scatter(ds1x, ds1y, color="grey")
plt.scatter(ds2x, ds2y, color="black")
plt.show()
