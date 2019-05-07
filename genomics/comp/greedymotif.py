"""
Greedy motif search scans each DNA sequence only once. Once we have scanned a 
particular sequence i, we decide which of its l-mer has the best contribution 
to the partial alignment score Score(s,i,DNA) for the first i sequences and 
immediately claim that this l-mer is part of the alignment.
"""

def gms(s, t, n, l):

