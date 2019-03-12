import numpy as np

"""
The dynamic programming approach by Michael S. Waterman, Temple F. Smith and William A. Beyer (1976)
computes optimal global alignments for two sequences and allows an arbitrary scoring of consecutive
gaps (insertions/deletions). Thus, it can incorporate the assumption that a single large insertion/deletion
event is biologically more likely to happen compared to many small insertions/deletions by using e.g. a
logarithmic or affine gap scoring function.


"""
