import re

"""
Greedy motif search scans each DNA sequence only once. Once we have scanned a 
particular sequence i, we decide which of its l-mer has the best contribution 
to the partial alignment score Score(s, i, DNA) for the first i sequences and 
immediately claim that this l-mer is part of the alignment.

For whatever scoring metric we use, the underlying motif search algorithm can
very based on different ways of searching and identifying for a given sequence
string.
"""

def regms(s):
    """
    Something similar to what is described above using regex.
    Other potential answers can be explored.
    Return repeated motifs in s.
    """
    m = re.search(r'(.+)\1+', s) 
    result = [] # output result
    prevind = 0 # previous index
    if m: # as long as there are repeats
        i, j = m.span() # m.start(group) and m.end(group)
        sub = s[i:j] # find the first subsequence
        ind = (sub + sub).find(sub, 1) # set the index
        sub = sub[:ind] # get the next one
        if len(sub) > 1:
            result.append([sub, (1 + prevind + 1, j + prevind + 1)])
        prevind += j # move the previous index up one
        return regms(s[j:], prevind) # run it again
    return result

def probgms(s):
    """
    Motif search using probability matrices.
    """
