import numpy as np

"""
Simple threshold clustering seeks to consider items in a way that's sensitive
to the arrangement of data points.
"""

def euclideanDist(vectorA, vectorB):
    diff = vectorA - vectorB
    
