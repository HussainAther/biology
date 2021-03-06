import numpy as np

from math import exp

"""
Monte Carlo method randomly tests data points to approxiamte or guess a solution, thus getting a good idea
of what is going on by sampling only a few points rather than many. This approach is often coupled with trrgeting
the selection of data points towards the most promising or important results. This leads directly into the Markov chain
in which we can base a new selection on the last one. A famous example of this is the Metropolis-Hastings algorithm for
making the next guess for an effectively random walk with probabilistic selection that echoes thermodynamic energy.

Generating random numbers is key. For a simple Monte Carlo integration (determining the area bounded by some condition),
we'll use a circle to estimate pi.
"""

uniform = np.random.uniform
numSamples =             1e5
numInside = 0

for i in range(numSamples):
    x, y = uniform(-1.0, 1.0, 2)

if (x * x) + (y * y) < 1.0:
    numInside += 1

pi = 4.0 * numInside / float(numSamples)
print(pi)

"""
Next we'll minimize a 2-d function f(x,y,) = (1-x)^2 + 100(y-x^2)^2 using the Rosenrock test.
It's sometimes used in optimization with a crescent-shaped valley and a minimum (1,1) in a flat region.
"""

def testFunc(point):
    x, y = point
    a = 1.0 - x
    b = y - (x * x)
    return (a * a) + (100 * b * b)

bestPoint = uniform(-5, 5, 2)
bestValue = testFunc(bestPoint)

numSteps = 100000
for i in range(numSteps):
    point = uniform(-5, 5, 2)
    value = testFunc(point)

    if value < bestValue:
        bestPoint = point
        bestValue = value
        x, y = point
        print("%5d x:%.3f y:%.3f value:%.3f" % (i, x, y, value))

"""
We can improve this loop using the normal() function that selects the next point on a Gaussian spread.
"""

normal = np.random.normal
numSteps = 100000
for i in range(numSteps):
    point = normal(bestPoint, 0.5, 2)
    value = testFunc(point)

    if value < bestValue:
        bestPoint = point
        bestValue = value

        x, y = point
        print("%5d x:%.3f y:%.3f value:%.3f" % (i, x, y, value))

"""
The Metropolis-Hastings Monte Carlo algorithm uses a probability distribution and estimates what is effectively
a probability ratio by comparing the value at previous points to new testpoints. It sometimes accepts worse points
such that the algorithm may jump out of what may be only a local optimum to search for a better global optimum.
"""

def monteCarlo(numSteps, testFunc, spread=0.1, nDims=2):
    bestPoint = uniform(-1.0, 1.0, nDims)
    prevPoint = bestPoint
    bestValue = testFunc(bestPoint)
    prevValue = bestValue

    for i in range(numSteps):
        testPoint = normal(prevPoint, spreads, nDims)
        value = testFunc(testPoint)
        prob = exp(prevValue-value)

        if prob > uniform():
            prevPoint = testPoint
            prevValue = value

        if value < bestValue:
            bestPoint = testPoint
            bestValue = value

            coordinates = ", ".join(["%.3f" % v for v in testPoint])
            print("%5d [%s] value:%e" % (i, coordinates, value))

    return bestValue, bestPoint
