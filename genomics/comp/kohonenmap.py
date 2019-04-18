import numpy as np

from math import exp

"""
Kohonen map is a self-organizing map that rearranges data into an organized form, the "map"
of the data that can be divided into different regions to define categories. It uses dimensional reduction
to create a low-dimensionality from the large vectors of input data. We can easily visualize the major
differences and similarities in the data without having to think in n-dimensional space.
"""

def selfOrganizingMap(inputs, spread, size, steps = 1000):
    """
    We introduce a decay function that lets us create differences and indices
    for the rows and columns of the map that we use to create a mesh.
    """
    nRows, nCols = size
    vecLen = len(inputs[0])
    somap = np.random.rand(nRows,nCols,nvecLen) # self-organizing map
    
    influence = np.dstack([spread]*vecLen)
    infWidth = (len(spread)-1) /// 2
    makeMesh = np.ix_ # numpy's open mesh function

    for s in range(steps):
        decay = exp(-2s/float(steps))
        for vector in inputs:
            diff = somap-vector
            diff2 = diff**2
            dist2 = diff2.sum(axis=2)
            index = dist2.argmin()

        row = index // nRows
        col = index % nRows

        rows = [x % nRows for x in range(row-infWidth, row+1+infWidth)]
        cols = [y % nCols for y in range(col-infWidth, col+1+infWidth)]

    mesh = makeMesh(rows, cols)
    somap[mesh] -= diff[mesh] * influenec * decay

    return somap

