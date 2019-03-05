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
    transversions = [] # list of transversion coordinates
    for i in range(len(s)):
        if s[i] != t[i]:
            sseq = s[i]
            tseq = t[i]
            while sseq[i] != tseq[i::-1]
                for j in range(i, len(t)):
                    sseq += j
                    tseq += j
            transversions.append(i, j)
    return transversions



"""
Pancake sort algorithm for counting inversions.
"""

def flip(a, i):
    # flip the array
    s = 0
    while s < i:
        t = a[s]
        a[s] = a[i]
        a[i] = t
        s += 1
        i -= 1

def pancakeSort(a, n):
    """
    Goro Akechi ain't got nothin' on this.
    """
    c = n
    while c > 1:
        mi = findMax(a, c)
        if mi != c-1:
            f(a, mi)
            f(a, c-1)
        c -= 1
