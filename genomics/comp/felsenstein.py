"""
Felstenstein's algorithm (felsenstein) for likelihood of phylogenetic trees.
"""

class Node:
    """
    Initialize each node in the tree.
    """
    def __init__(self, position, leaf):
        self.position = position
        self.isLeaf = leaf 
