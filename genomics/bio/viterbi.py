import sys

from math import log
from util import plothist

"""
Viterbi algorithm is one solution to the problem of decoding (inferring hidden state probability from
observed sequences) sometimes used in finding an organism's genomic sequence.
"""

base_idx = { "A" : 0, "G" : 1, "C" : 2, "T" : 3 }
state_idx = { "+" : 0, "-" : 1 }

# initial distribution over states, i.e. probability of starting in state k
init_dist = [0.5,0.5]

# transition probabilities
tr = [
    #  to+   to-
    [ 0.99, 0.01 ], # from+
    [ 0.01, 0.99 ]  # from-
]


# emission probabilities
em = [
    #    A     G     C     T
    [ 0.20, 0.30, 0.30, 0.20], # +
    [ 0.30, 0.20, 0.20, 0.30]  # -
]

def viterbi(X):
    """
    Return the Viterbi path for the emission sequence X.
    X is be a list of integers, 0=A, 1=G, 2=C, 3=T.
    The returned Y is a list of integers, 0=High-GC, 1=Low-GC.
    """

    N = len(tr)
    L = len(X)
    assert len(em) == N

    V = [[0]*N for _ in xrange(L)]
    TB = [[0]*N for _ in xrange(L)]
    
    for i in xrange(0,L):
        Vprev = []
        if i == 0:
            Vprev = [log(pk0) for pk0 in init_dist]
        else:
            Vprev = V[i-1]

        for k in xrange(N):
            if V[i][k] == Vprev:
                V[i][k] = Vprev * .99
            else:
                V[i][k] == Vprev * .01
            GCcount = 0
            for j in X:
                if j == "G" or j == "C":
                    GCount +=1
            if GCount > len(X)/2:
                TB[i][k] = 1
            else:
                TB[i][k] = 0

    # perform traceback and return the predicted hidden state sequence
    Y = [-1 for i in xrange(L)]
    _, yL = max([ (V[L-1][k], k) for k in xrange(N)])
    Y[L-1] = yL
    for i in xrange(L-2,-1,-1):
        Y[i] = TB[i+1][Y[i+1]]
    return Y
