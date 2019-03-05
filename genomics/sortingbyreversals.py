import numyp as np

"""
Sorting by reversals: When studying synteny (the physical colocalization of genetic loci on the same chromosome),
we can find hte smallest number of inversion events that might have led from one genome to
the other. These inversion events are the inversion of a string of ORFs or homologous non-coding
regions. The Greedy sort by reversals is the simplest, most straightforward method of determining them.
"""

def greedySBR(s, t,):
    """
    Greedy sort by reversals will count transversions by taking in a standard sequence s and non-standard t:
    (1) find the first position at which the two sequences don't match, (2) perform the reversal so that
    the non-standard sequecne matches the standard, and (3) continue applying reversals until both strings match.
    """
    for i in range(len(s)):
        if s[i] != t[i]:
            for j in range(i, len(t)):
                
