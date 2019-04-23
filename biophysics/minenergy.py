"""
RNA has some of the properties of both DNA and proteins. It has the same 
information storage capability as DNA due to its sequence of nucleotides. 
But its ability to form 3D structures allows it to have enzymatic properties 
like those of proteins. We call the linear polynucleotide chain the primary 
structure of RNA; the intra-molecular base-pairing is the secondary structure.

Secondary structures are important for RNA functions. For example, tRNA 
molecules all form similar structures (a characteristic "clover-leaf", 
one branch of which is the "anticodon loop", which pairs with a codon as 
part of the process of translating a messenger RNA into protein). Overall, 
the structures are critical to their function.

Of all possible 16 base-pairings, 6 of them are considered stable. The Watson-Crick 
base pairs, C-G and A-U, form stable base pairs with each other through 
the creation of hydrogen bonds between donor and acceptor sites on the bases. 
In addition, the G-U wobble pair is weaker but still stable. Other bases also 
pair sometimes, especially if chemically modified.

We can find the strucure with minmal free energy. We can look at hairpin, bulge,
interior, stacking, and multi loops (multi-loop).
"""

def stack(x, y):
    """
    Return stacking energy between two base pairs x and y.
    """
    return abs(x-y)

def asymmetry(x, y):
    """
    Asymmetry function gives a penalty for asymmetric loops.
    """
    return -abs(x-y)

def eL(a, b, c, d):
    """
    Energy of bulge or internal loop with exterior base pair
    (a, b) and interior base pair (c, d).
    """
    return len(c-a+d-b-2) + stack(a, b) + stack(c, d) + asymmetry(c-a-1, b-d-1) 
    
def eM(args):
    """
    For a list of tuples args, determine the energy of an 
    optimal structure of the subsequence from i through j 
    in which (i, j) closes a multibranched loop.
    """ 
    i0 = args[0][0]
    i1 = args[1][0]
    j0 = args[0][1]
    jl = args[-1][1]
    summa = 0
    for k in range((len(args)//2) - 1):
        summa += (args[k+1][0] - args[k][1] -1)
    return i1-i0-1+j0-jl+summa      
 
def VBI(i, j, l): 
    """
    
    """  
    return min(e
