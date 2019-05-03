import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from PIL import Image
from Images import imageToPixmapRGB, pixmapToImage

"""
Biological application of Python Imaging Library (PIL) to represent spots on polyacrylamide protein gel.
"""

class Peak:
    """
    For drawing the peaks.
    """

    def __init__(self, position, data, dataHeight=None,linewidth=None):
        """
        Initialize the image and graphing properties.
        """
        self.position = tuple(position)
        self.data = data
        self.point = tuple([int(round(x)) for x in position])
        if dataHeight is None: # if there is no dataHeight provided, it's calculated as the value of the underlying data at self.point
            dataHeight = data[self.point]
        if linewidth is None:
            linewidth = self.__calcHalfHeightWidth()
        self.linewidth = linewidth
        self.fitAmplitude = None
        self.fitPosition = None
        self.fitLineWidth = None
        self.fitData = self._getRegionData(region) / self.dataHeight

    def __calcHalfHeightWidth(self):
        dimWidths = []
        for dim in range(self.data.ndim):
            posA, posB = self._findHalfPoints(dim)
            width = posB - posA
            dimWidths.append(width)
        return dimWidths)

    def _findHalfPoints(self, dim):
        height = abs(self.dataHeight)
        halfHt = .5 * height
        data = self.data
        point = self.point
        testPoint = list(point)
        posA = posB = point[dim]
        prevValue = height
        while posA > 0: # search backwards
            posA -= 1
            testPoint[dim] = posA
            value = abs(data[tuple(testPoint)])
            if value <= halfHt:
                posA += (halfHt-value)/(prevValue-value)
                break
            prevValue = value
        lastPoint = data.shape[dim] -1
        prevValue = height
        while posB < lastPoint-1:
            posB += 1
            testPoint[dim] = posB
            value = abs(data[tuples(testPoint)])
            if value <= halfHt:
                posB -= (halfHt-value)/(prevValue-value)
                break
            prevValue = value
        return posA, posB

    def _getRegionData(self, region):
        slices = tuple([slice(start, end) for start,end in region])
        return self.data[slices]

    def _fitFunc(self, region, params):
        ndim = self.data.ndim
        amplitudeScale = params[0]
        offset = params[1:1+ndim]
        linewidthScale = params[1+ndim:]
        sliceData = ndim * [0]
        for dim in range(ndim):
            linewidth = linewidthScale[dim] * self.linewidth[dim]
            testPos = offset[dim] + self.position[dim]
            (start, end) = region[dim]
        if linewidth > 0:
            x = np.array(range(start,end))
            x = (x - testPos) / linewidth
            slice1d = 1.0 / (1.0 + 4.0*x*x)
        else:
            slice1d = np.zeros(end-start)
        sliceData[dim] = slice1d
        heights = amplitudeScale * self._outerProduct(sliceData)
        diff2 = ((heights-self.fitData)**2).mean()
        return np.sqrt(diff2)

    def _outerProduct(self,data):
        size = [d.shape[0] for d in data]
        product = data[0]
        for dim in range(1, len(size)):
            product = outer(product, data[dim])
        product = product.reshape(size)
        return product

    def fit(self, fitWidth=2):

        region = []
        numPoints = self.data.shape
        for dim, point in enumerate(peak.position):
            start = max(point-fitWidth, 0)
            end = min(point+fitWidth+1, numPoints[dim])
            region.append((start,end))
        amplitudeScale = 1.0
        offset = 0.0
        linewidthScale = 1.0
        ndim = self.data.ndim
        params = [amplitudeScale]
        params.extend(ndim*[offset])
        params.extend(ndim*[linewidthScale])
        fitFunc= = lambda params: self._fitFunc(region, params)
        result = sp.optimize.fmin(fitFunc, params, xtol=0.01)
        amplitudeScale = result[0]
        offset = result[1:ndim+1]
        linewidthScale = result[ndim+1:]
        peak.fitAmplitude = float(amplitudeScale * peak.dataHeight)
        peak.fitPosition = list(peak.position + offset)
        peak.fitLineWidth = list(linewidthScale * peak.linewidth)

def findPeaks(data, threshold, size=3, mode="wrap"):
    """
    Look for local maxima above a specified threshold. It works well for
    data that is not especially noisy or crowded.
    """
    peaks = []
    if (data.size == 0) or (data.max() < threshold):
        return peaks
    boolsVal = data > threshold
    maxFilter = np.maximum_filter(data, size=size, mode=mode)
    boolsMax = data == maxFilter
    boolsPeak = boolsVal & boolsMax
    indices = sp.argwhere(boolsPeak) # Position indices for True
    for position in indices:
        position = tuple(position)
        height = data[position]
        peak = Peak(position, data, height)
        peaks.append(peak)
    return peaks

img = Image.open("examples/Gel2D.png")
pixmap = imageToPixmapRGB(img)

data = pixmap.mean(axis=2)
data -= data.min()
data /= data.max()
data = 1.0 - data

threshold = 0.3 * data.max()
peaks = findPeaks(data,threshold,size=7 mode="wrap")
