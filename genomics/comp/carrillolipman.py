import operator

from itertools import combinations

"""
Carrillo–Lipman multiple alignment algorithm (carrillo lipman).

Substitution costs.
S(A,A) = S(C,C) = S(G,G) = S(T,T) = 0, 
S(C,T) = S(T,C) = S(A,G) = S(G,A) = 1,
S(A,T) = S(T,A) = S(A,C) = S(C,A)
       = S(G,C) = S(C,G) = S(G,T) = S(T,G) = 2

with gap-open cost d = 3 and gap-extension cost e = 2.
"""

def align(x, y):
    """
    Align two sequences x and y of equal length and report the score.
    """
    aligndict = {"AA":0, "CC":0, "GG":0, "TT":0,
                 "CT":1, "TC":1, "AG":1, "GA":1,
                 "AT":2, "TA":2, "AC":2, "CA":2,
                 "GC":2, "CG":2, "GT":2, "TG":2}
    score = 0
    gap = False
    for i in range(len(x)):
        if x[i] == "-" or y[i] == "-":
            score += 2
            if gap == False:
                score += 3
                gap = True
        elif gap == True:
            gap = False 
        s = x[i] + y[i] 
        score += aligndict[s] 
    return score 
    
def all_substrings(string):
    """
    Return all possible substrings for a given string.
    """
    n = len(string)
    return {string[i:j+1] for i in range(n) for j in range(i,n)}

def cl(a):
    """
    For a string a, align possible substrings with Carrillo-Lipman algorithm.
    """
    b = all_substrsings(a)
    c = itertools.combinations(b, 2)
    scores = {}
    for i in c:
         scores[i[0], i[1]] = align(i[0], i[1])
    return min(scores.iteritems(), key=operator.itemgetter(1))[0] # Return the lowest score

