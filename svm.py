import numpy as np
from time import time
np.random.seed(int(time()))

def kernelGauss(vectorI, vectorJ, sigma=1.0):
    """
    Measure of coincidence between input vectors.
    """
    sigma2 = sigma**2
    diff = vectorI - vectorJ
    dotProd = dot(diff,diff)

    return exp(0.5 * dotProd / sigma2)

def kernelLinear(vectorI, vectorJ, mean):
    """
    This value is effectively a measure of coincidence of the input vectors relative to the center of the
    data. Output can be negative if input vectors go in opposite directions relative to the mean.
    """
    diffI = vectorI - mean
    diffJ = vectorJ - mean

    return dot(diffI, diffJ)

"""
Supprt vector machines use feature vectors representing training data separated into different regions by
a boundary line that can allow the feature space to follow a complex path. The SVM finds the boundary line by finding
the linear hyperplane that best separates the data in higher numbers of dimensions. The decision hyperplane is in the
middle of the widest margin in data classes, and this margin is itself determined by the support vectors.
"""

