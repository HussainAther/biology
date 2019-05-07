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

def probgms(a, k):
    """
    Motif search using probability matrices. For array of strings a,
    use a profile matrix from the first k-mer. The probability matrices
    let us use the best probabilities from the profile matrices of the best
    k-mer groups in each case.  
    """
    best_motifs = [i[0:k] for i in a] # Initialize list of best motifs
    score = 0 
    c = 0
   


def matrixgen(s, matrix):
    """
    Generate a DNA matrix for a string s if there isn't one already for the input 
    matrix. If there is one, sum up the probability values and average them
    to make sure they're appropriate. 
    """
    if len(matrix)==0:
        for ch in s: # for each character in the string
            new_dict = {"A":0, "C":0, "G":0, "T":0}
            new_dict[ch] = 1
            matrix += [new_dict]
        print(matrix)
        return matrix
    elif len(matrix) == 3:
        for j in range(len(s)):
            matrix[j][s[j]] += 1
        Sum = [sum(matrix[0].values()), sum(matrix[1].values()), sum(matrix[2].values())]
        for columns in range(len(matrix)):
            for keys in matrix[columns]:
                matrix[columns][keys] /= Sum[columns]
        print(matrix)
        return matrix  
