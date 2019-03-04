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
