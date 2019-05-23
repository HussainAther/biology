"""
Felstenstein's algorithm (felsenstein) for likelihood of phylogenetic trees.
"""

class Node:
    """
    Initialize each node in the tree.
    """
    def __init__(self, position, leaf):
        """
        Initialize each position in the tree
        with a boolean of whether it is a leaf.
        """
        self.position = position
        self.isLeaf = leaf
        
def felsenstein(a, x):
    """
    For an input dictionary a of position nodes and and which nodes they connect to and
    x array of residues, return the likelihood sites for the sequences.
    """ 
    n = len(a) # length of dictionary a
    k = 2*n -1 # possible node children
    p = [] # probabilities
    while k > 0:
        if a[k].isLeaf: # if k is a leaf
            if a[k] == x[k]:
                p.append(1)
            else:
                p.append(0)
        else: # if k is not a leaf we can speculate on two daughter nodes that may follow from it
            p.append(p[k+1] * p[k+2])
    return p
