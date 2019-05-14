import sys
import fileinput
import math

"""
Given a set of putative exons, find a maximum set of nonoverlapping putative exons.
We use dynamic programming to solve the Exon Chaining (exon chaining) problem across various 
intervals. For a graph G with 2n vertices, n which represents the starting interval position
and the ending position. For the intervals in increasing order, we get an array of vertices

"""

import sys
import fileinput
import math

class Node:
    """
    Define nodes on the graph.
    """
    def __init__(self, location, score, startNode = None):
        self.location = location
        self.score = score
        self.startNode = startNode
        self.isEnd = self.startNode != None

def buildGraphFromInput(lines):
    """
    For a list of lines l, we initialize the node
    using the start and end to construct a graph.
    """
    graph = []
    #read file and initializes nodes in graph
    for line in lines:
        if len(line.strip()) == 0:
            continue

        #parse line
        start,end,score = [int(num.strip()) for num in line.split()]

        #initialize node
        startNode = Node(start, score)
        endNode = Node(end, score, startNode)

        #add nodes to graph
        graph.extend([startNode, endNode])

    #sort graph and initialize index references on each node
    graph.sort(key=lambda a: (a.location, a.startNode.location if a.startNode else -1))

    for i,n in enumerate(graph):
        n.index = i

    return graph

def getScore(graph):
    scores = [0]

    for node in graph:
        if not node.isEnd:
            scores.append(scores[-1])
        else:
            scoreIndex = node.startNode.index
            score = node.score

            newScore = max(scores[scoreIndex] + score, scores[-1])
            scores.append(newScore)

    return scores[-1]
