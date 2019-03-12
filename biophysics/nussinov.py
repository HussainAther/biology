import sys
import numpy as np

"""
We implement that Nussinov RNA folding algorithm for predicting the secondary structure of RNA.
We get a two-dimensional secondary structure of RNA by determining how pairing between nucleotides
in the primary sequence give rise to a two-dimensional architecture.

We must obtain an optimal pairing for a string of nucleotides. In order to know how good our
structure is, we assign a score to it. One possible scoring scheme could be adding 1 to the score per paired
set of nucleotides, and 0 otherwise. So in other words, we want a pairing that will give us the highest possible
score.

We can write this quantity as optimal(i,j) where i and j are the indices of the sequence between which we
obtain the pairing score. Our algorithm is therefore going to compute a folding score for all substrings
bound by i and j and store the value in what is known as a dynamic programming table.

Our dynamic programming table will be a N x N array where N is the length of our sequence. So now that we
have a way of measuring how good a structure is, we need a way to evaluate scores given a subsequence.

To do this, we set some rules on the structure of an RNA sequence:

If i and j form a pair:

1. The pair i and j must form a valid watson-crick pair.
2. i<j−4. This ensures that bonding is not happening between positions that are too close to each other, which would produce steric clashes.
3. If pair (i,j) and (k,l) are in the structure, then i<k<j<l. This ensures that there is no crossing over of pairs which would result in pseudoknots.
4. No base appears in more than one pair.
"""

def pair(t):
    """
    Check if the pair of nucleotides bond with one another.
    """
    if t in [("A", "U"), ("U", "A"), ("C", "G"), ("G", "C")]:
        return True
    return False

def optimal(i, j, seq):
    """
    Return optimal pairing score between two indices i and j.
    """
    if i >= j - 4: # no pairs when i and j are fewer than four bases apart
        return 0
    else:
        u = optimal(i, j-1, sequence) # recursively look for optimal pairing score among unpaired possibilities
        p = [1 + optimal(i, t-1, seq) + optimal(t+1, j-1, sequence) for t in range(i, j-4) if pair((seq[t], seq[j]))] # among paired options

    if not p:
        p = [0]
    p = max(p)

    return max(u, p)
