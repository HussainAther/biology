import numpy as np

from random import random

"""
The Negative Selection algorithm is based off the self-nonself discrimintaino behavior
in the mammalian acquired immune system. The clonal selection theory of acquired immunity
accounts for the adaptive behavior of the immune system that includes the ongoing selection
and proliferation of cells that select for potentially harmful material in the body. We use 
negative selection to select for self-reactive cells during cell creation and cell proliferation.
We observe this in the creatino of T-lymphocytes in the thymus.
"""

def randvec(minmax):
    """
    Return a random vector from minmax.
    """
    result = []
    for i in minmax:
        result.append(minmax[i][0] + ((minmax[i][1] - minmax[i][0]) * random())) 
    return result

def eucdist(a, b):
    """
    Return Euclidean distance between two vectors a and b.
    """
    return np.linalg.norm(a-b) 
 
def contains(vec, space):
    """
    Return True if the vector vec is in space. False otherwise.
    """
    for i in vec:
         for j in i:
               if i< space[j][0] or i > space[j][1]:
                   return False
    return True

def matches(vec, dataset, min_dist):
    """
    If there are any matches between the data set and the vector
    based off some minimum distance.
    """
    for p in dataset:
        dist = eucdist(vec, p[:vec])
            if dist <= min_dist:
                 return True
    return False

def gendet(maxdet, searchspace, selfdataset, mindist):
    """
    Return the detectors for the search space, dataset, and minimum distance between.
    """
    det = []
    while len(det) < maxdet:
        vec = randvec(searchspace) 
        det = det[:vec]
        if not matches[det[:vec], selfdataset, mindist):
            if not matches(det[:vec], det, 0):
                det.append(det) 
    return det


def genselfdata(num, selfspace, searchspace):
    """
    Generate the dataset selfdata for our algorithm.
    """
    selfdata = []
    while len(selfdata) < num:
        p = randvec(searchspace)
        if matches(p, selfdata, 0) and contains(p, selfspace):
            selfdata.append(p)
    return selfdata

def applydet(det, bounds, dataset, mindist, trials=50):
    """
    Apply our detectors to our algorithm.
    """
    correct = 0
    for i in range(trials):
        input = randvec(bounds)
        actual = []
        expected = []
        correct = 0
        for j in input:
            actual.append(matches(input, det, mindist)
            expected.append(matches(input, dataset, mindist)
            if actual[-1] == expected[-1]: # if actual and expected values are the same
                correct += 1
    return corret

def runit(bounds, selfspace, maxdetect, maxself, mindist):
    """
    Execute the algorithm.
    """
    selfdataset = genselfdata(maxself, selfspace, bounds)
    detect = gendet(maxdetect, bounds, selfdataset, mindist)
    applydet(detect, bounds, selfdataset, mindist)
    return detect
