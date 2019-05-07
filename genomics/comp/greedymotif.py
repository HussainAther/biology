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
    m = re.search(r'(.+)\1+', s)
    result = [] # output result
    prevind = 0 # previous index
    if m:
        i, j = m.span() # m.start(group) and m.end(group)
        sub = s[i:j]
        ind = (sub + sub).find(sub, 1)
        sub = sub[:ind]
        if len(sub) > 1:
            result.append([sub, (1 + prevind + 1, j + prevind + 1)])
        prevind += j
        return gms(s[j:], prevind) 
    return result
