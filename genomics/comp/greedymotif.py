import re

"""
Greedy motif search scans each DNA sequence only once. Once we have scanned a 
particular sequence i, we decide which of its l-mer has the best contribution 
to the partial alignment score Score(s, i, DNA) for the first i sequences and 
immediately claim that this l-mer is part of the alignment.
"""

def gms(s):
    """
    Something similar to what is described above using regex.
    Other potential answers can be explored.
    Return repreated motifs in s.
    """
