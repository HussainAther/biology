import numpy as np

"""
Extract innate trends within data. The principal components (pca) are the vectors that give the
best separation of data items in terms of their covariances.
"""

def principalComponentAnalysis(data, n=2):
    """
    Steps of principal component analysis
    """
    samples, features = data.shape
    meanVec = data.mean(axis=0)
    dataC = (data - meanVec).T # take the tranpose of the data with the mean vectors
    covar = np.cov(dataC) # covariance
    evals, evecs = np.linalg.eig(covar) # eigenvalues and eigenvectors
    indices = evals.argsort()[::-1]
    evecs = evecs[:,indices]
    basis = evecs[:,:n]
    energy = evals[:n].sum()
    return basis, energy

def extractPrincnipalComponent(data, precision=1e-9):
    """
    Extract the principal components
    """
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
    """
    Linear discriminant analysis
    """
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
