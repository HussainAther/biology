import numpy as np

"""
The dynamic programming approach by Michael S. Waterman, Temple F. Smith and William A. Beyer (1976)
computes optimal global alignments for two sequences and allows an arbitrary scoring of consecutive
gaps (insertions/deletions). Thus, it can incorporate the assumption that a single large insertion/deletion
event is biologically more likely to happen compared to many small insertions/deletions by using e.g. a
logarithmic or affine gap scoring function.

To this end, all possible gap lengths are explicitly considered within the recursion, which increases the time complexity from
O(n^2) to O(n^3) when compared to the approach by Needleman and Wunsch (1970) that applies a linear gap scoring. Under the assumption that both input sequences
a and b stem from the same origin, a global alignment tries to identify matching parts and the changes needed to transfer one sequence into the other.
The changes are scored and an optimal set of changes is identified, which defines an alignment. The dynamic programming approach
tabularizes optimal subsolutions in matrix D. An entry D_ij represents the best score for the alignment of the prefixes a_1i with b_ij. We use adding recursion.


"""
