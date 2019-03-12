import numpy as np
import sys

"""
In 1990 John S. McCaskill introduced an efficient dynamic programming algorithm to compute
the partition function (Z=∑Pexp(−E(P)/RT)) over all possible nested structures P that can be
formed by a given RNA sequence S with E(P) = energy of structure P, R = gas constant, and
T = temperature.
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
