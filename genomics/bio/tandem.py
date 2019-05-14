"""
A tandem repeat Pk of a pattern P = p1 . . . pm is a pattern of length k ·m formed
by concatenating k copies of P. The tandem repeat problem is to find the best
global alignment between some interval of a text T = t1 . . . tn and a tandem repeat
Pk for some k. When computing global alignments the premium for matching is
+1 and the mismatch and gap penalties are −1.

For example given the pattern P = GGT and the text T = ACCGGTGCTGTAA the best
tandem repeat is given by k = 3 and the following alignment:

ACCGGTGCTG-TAA
   GGTGGTGGT
"""

def tandem(p, t, score):
    """
    Fill an array D with normal recursion for the tandem repeat alignment problem.
    For some pattern p and text t, we determine the best tandem repeat when aligning
    p to t.
    """ 
    for i in p:
        if i < p: 
            score -= 1
