"""
RNAs adopt complex three-dimensional structures that are important for many biological functions. 
Pairs of positions in RNA with complementary nucleotides can form bonds. Bonds (i, j) and (i′, j′) 
are interleaving if i < i′ < j < j′ and noninterleaving otherwise (fig. 6.30). Every set of 
noninterleaving bonds corresponds to a potential RNA structure. In a very naive formulation of the 
RNA folding problem, one tries to find a maximum set of noninterleaving bonds. The more adequate model, 
attempting to find a fold with the minimum energy, is much more difficult.
"""

def noninter(a, b):
    """
    For two lists of RNA positions a and b, find the maximum set of noninterleaving bonds.
    """
    bonds = []

