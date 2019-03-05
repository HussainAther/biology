import numpy as np

"""
Extract innate trends within data. The principal components are the vectors that give the
best separation fo data items in terms of their covariances.
"""

def principalComponentAnalysis(data, n=2):
    samples, features = data.shape

    meanVec = data.mean(axis=0)
    dataC = (data - meanVec).T

    covar = np.cov(dataC)
    evals, evecs = np.linalg.eig(covar)

    indices = evals.argsort()[::-1]
    evecs = evecs[:,indices]

    basis = evecs[:,:n]
    energy = evals[:n].sum()

    return basis, energy

def extractPrincnipalComponent(data, precision=1e-9):
    samples, features = data.shape
    meanVec = data.mean(axis=0)
    dataC = data - meanVec

    pc1 = random.random(features)
    pc0 = pc1 - 1

    while abs((pc0-pc1).sum()) > precision:
        t = np.zeros(features)
        for datum in dataC:
            t += np.dot(datum, pc1) * datum

        pc0 = pc1
        pc1 = t / np.sqrt(np.dot(t,t,))

    return pc1

"""
Linear discriminant analysis finds the matrix that maximises the scatter (separation) between
data sets relative to the separation within each data set. The separation within is the weighted
covariances of teh data sets separately and the separation betwen is the difference between
their means.
"""

def twoClassLda(dataA, dataB):
    meanA = dataA.mean(axis=0)
    meanB = dataB.mean(axis=0)

    covA = np.cov(dataA.T)
    covB = np.cov(dataB.T)

    nA = len(dataA) - 1
    nB = len(dataB) - 1

    scatterWithin = nA * covA + nB * covB
    scatterBetween

    discrim = np.dot(np.linalg.inv(scatterWithin), scatterBetween)

    transfA = np.dot(dataA, discrim.T)
    transfB = np.dot(dataB, discrim.T)

    divide = np.dot(discrim, (meanA + meanB))/2

    return transfA, transfB, divide
