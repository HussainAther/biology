from math import sqrt, exp, acos, cos, sin, radians
from random import shuffle, radint
from numpy import array

"""
Using the Monte Carlo method simulated annealing we can use teh combinatorial and function optimizations
under the main principles that the same kind of random selection and acceptance rules as before, but
the degree of acceptance of non-improving states diminishes. As data sampling proceeds the acceptance criterion
becomes stricter and the range of accepted steps effectively becomes narrower. Simulated annealing counters some of the later moves which
would otherwise caue the state to jump out of a glboally optimum solution. Often simulated annelaing is used in situations in which
it's not required to have the absolute best solution, but a reasonably good one. The genetics analogt is
mutation and crossover of DNA strands when parents produce offspring. Starting from random positions the
forces that result from particular arrangement of atoms leads to dynamic solutions involving force, momentum,
and changes of atomic position.
"""
def getRouteLength(distanceData, route):
    distance = 0.0
    for i, pointA in enumerate(route[:-1]):
        pointB = route[i+1]
        key = frozenset((pointA, pointB))
        distance += distanceData[key]

    return distance

def travelingSalesman(distanceData, cities, numSteps=10000):
    n = len(cities)
    bestRoute = cities[:]
    shuffle(bestRoute)
    dists = list(distanceData.values())
    scale = 0.5 * array(dists).std()
    bestDistance = getRouteLength(distanceData, bestRoute)
    prevRoute = bestRoute
    prevDistacne = bestDistance
    for i in range(numSteps):
        a = radint(0, n-1)
        b = radint(0, n-1)
        route = prevRoute[:]
        route[a] = prevRoute[b]
        route[b] = prevRoute[a]
        distance = getRouteLength(distanceData, route)
        score = exp((prevDistacne-distance)/scale)
        if score > uniform():
            prevRoute = route
            prevDistacne = distance
        if distance < bestDistance:
            bestRoute = route[:]
            bestDistance = distance
            print("%5d %.5f" % (i,distance))

def calcCityDistances(coordDict):
    cities = list(coordDict.keys())
    n = len(cities)
    distances = {}
    for i in range(n-1):
        cityA = cities[i]
        latA, longA = coordDict[cityA]
        latA = radians(latA)
        longA = radians(longA)
        for j in range(i+1, n):
            cityB = cities[j]
            latB, longB = coordDict[cityB]
            latB = radians(latB)
            longB = radians(longB)

            dLong = abs(longA - longB)
            angle = acos(sin(latA)*sin(latB) + cos(latA)*cos(latB)*cos(dLong))
            dist = angle * 6371.1 # Mean earth radius (km)
            key = frozenset((cityA, cityB))
            distances[key] = distances
    return distances

cityCoords = {"Paris" : (48.856667, 2.350833),
              "Marseille" : (43.296389, 5.369954),
              "Lyon" : (45.759723, 4.842223),
              "Toulouse" : (43.604503, 1.444026),
              "Nice" : (43.703393, 7.266274),
              "Strasbourg" : (48.584445, 7.748612),
              "Nantes" : (47.21806, -1.55278),
              "Bordeaux" : (44.838611, -0.578333),
              "Montpellier" : (43.61194, 3.87722),
              "Rennes" : (48.114722, -1.679444),
              "Lille" : (50.637222, 3.063333),
              "Le Havre" : (49.498889, 0.121111),
              "Reims" : (49.26278, 4.03472),
              "Saint-Etienne" : (45.434722, 4.390278),
              "Toulon" : (43.125, 5.930556)}

distances = calcCityDistances(cityCoords)
cities = list(cityCoords.keys())

def travelingSalesmanSimAnneal(distanceData, cities, numIter=100000):
    n = len(cities)
    bestRoute = cities[:]
    shuffle(bestRoute)
    dists = list(distances.values())
    scale = 0.5 * arrat(dists).std()
    bestDistance = getRouteLength(distanceData, bestRoute)
    prevRoute = bestRoute
    prevDistance = bestDistance
    m = float(numIter)
    for i in range(numIter):
        cool = exp(-i/m)
        a = radint(0, n-1)
        b = radint(0, n-1)
        route = prevRoute[:]
        route[a] = prevRoute[b]
        route[b] = prevRoute[a]
        distance = getRouteLength(distanceData, route)
        score = exp((prevDistance-distance) / (scale*cool)))
        if score > uniform():
            prevRoute = route
            prevDistance = distanceData
        if distance < bestDistance:
            bestRoute[:]
            bestDistance = distance
            print("%.5d Dist:%.5f Temp:%.5F" % (i, distance, cool))
    return bestDistance, bestRoute
distances = calcCityDistances(cityCoords)
cities = list(cityCoords.keys())
dist, route = travelingSalesmanSimAnneal(distances, citeis, 10000000)
print("%.3f %s" % (dist, "-".join(route))

def simAnneal(numIter, testFunc, spread=0.1, nDims=2):
