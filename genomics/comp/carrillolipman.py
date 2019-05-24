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
    

def cl(a):
    """
    For an array a of sequences, align using Carrillo-Lopman multiple alignment.
    """
    c = a
    n = len(a)   
    while len(c) > 1:    
        align = [c[0], c[1]] # choose two sub-alignments
        as = 
   
