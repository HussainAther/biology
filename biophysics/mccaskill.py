import numpy as np
import sys

"""
In 1990 John S. McCaskill introduced an efficient dynamic programming algorithm to compute
the partition function (Z=∑Pexp(−E(P)/RT)) over all possible nested structures P that can be
formed by a given RNA sequence S with E(P) = energy of structure P, R = gas constant, and
T = temperature.

Here, we provide a simplified version of the approach using a Nussinov-like energy scoring scheme,
i.e. each base pair of a structure contributes a fixed energy term Ebp independent of its context.
Given this, we populate two dynamic programming tables Q and Qbp. Qi,j provides the partition function
for subsequence from position i to j, while Qbpi,j holds the partition function of the subsequence
given the constraint that position i and j form a base pair (or 0 if no base pair is possible).
Watson-Crick as well as GU base pairs are considered complementary. The overall partition function
is given by Z=Q1,n for a sequence of length n.
"""

def mc():
    """
    Steps of the algorithm.
    """
    for k in range(1, L+1):
        for i in range(L-k+1):
            j = i + k -1

            sum = 0.0
            for h in range(i,j-1):
                sum += Q[i][h-1] * Eps_ij(h,j)


            Q[i][j] = Q[i][j-1] + sum

def init(N, q = 0.):
    """
    Initialize the dot bracket
    """
    dp = np.empty((N,N))
    dp[:] = np.NAN
    for i in range(N):
        dp[i][i] = q
        if i != 0:
            dp[i][i-1] = q
    return dp
