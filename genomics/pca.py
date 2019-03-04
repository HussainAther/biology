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

            
