from math import sqrt, exp
from random import shuffle, radint
from numpy import array

"""
How do you determine the best route between 15 of the largest cities in France?

Biological applications include its use in nuclear magnetic resonance as a common operation
to assign signals that occur in a spectrum to the amino acid residues of a protein chain.
"""

def getRouteLength(distanceData, route):
    distance = 0.0
    for i, pointA in enumerate(route[:-1]):
        pointB = route[i+1]
        key = frozenset((pointA, pointB))
        distance += distanceData[key]

    return distance
