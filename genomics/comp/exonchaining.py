import sys
import fileinput
import math

"""
Given a set of putative exons, find a maximum set of nonoverlapping putative exons.
We use dynamic programming to solve the Exon Chaining (exon chaining) problem across various 
intervals. For a graph G with 2n vertices, n which represents the starting interval position
and the ending position. For the intervals in increasing order, we get an array of vertices
in the graph. We have 3n - 1 edges (an edge between each li and ri fo weight wi for i from 
1 to n and 2n - 1 more edges of weight 0 that connect adjacent vertices (vi, vi+1) that all
make a path in the graph from v1 to v2n.
"""

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
    For a list of lines l of a garph, we initialize the node
    using the start and end to construct a graph.
    """
    graph = []
    # Read file and initializes nodes in graph
    for line in lines:
        # Parse line
        start, end, score = [int(num.strip()) for num in line.split()]
        # Initialize node
        startNode = Node(start, score)
        endNode = Node(end, score, startNode)
        # Add nodes to graph
        graph.extend([startNode, endNode])
    # Sort the graph and initialize index references on each node
    graph.sort(key=lambda a: (a.location, a.startNode.location if a.startNode else -1))
    for i, n in enumerate(graph):
        n.index = i
    return graph

def getScore(graph):
    """
    We get the score of the length of the longest path in the grpah
    that ends at the vertices vi.
    """
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
