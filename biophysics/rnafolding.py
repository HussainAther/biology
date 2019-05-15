"""
RNAs adopt complex three-dimensional structures that are important for many biological functions. 
Pairs of positions in RNA with complementary nucleotides can form bonds. Bonds (i, j) and (i′, j′) 
are interleaving if i < i′ < j < j′ and noninterleaving otherwise (fig. 6.30). Every set of 
noninterleaving bonds corresponds to a potential RNA structure. In a very naive formulation of the 
RNA folding problem, one tries to find a maximum set of noninterleaving bonds. The more adequate model, 
attempting to find a fold with the minimum energy, is much more difficult.
"""

def isNoninter(l):
    """
    For a list of two two-length tuple bond positions, determine if they're noninterleaving.
    """
    if l[0][0] < l[0][1] < l[1][0] < l[1][1]:
        return False
    elif l[1][0] < l[1][1] < l[0][0] < l[0][1]:
        return False
    return True 

def noninter(a, b):
    """
    For two lists of RNA positions a and b, find the maximum set of noninterleaving bonds.
    """
    bonds = []
    samplebond = []
    for i in a:
        for j in b:
            if samplebond == []:
                samplebond.append((i, j))
            elif len(samplebond) == 1:
                samplebond.append((i, j))
                if isNoninter(samplebond):
                    bond.append(samplebond)
                samplebond = []
    return samplebond   
   
