from math import sqrt, exp, acos, cos, sin, radians
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

dist, route = travelingSalesman(distances, cities, 1000000)
print("%.3f %s" % (dist, ", ".joint(route)))
