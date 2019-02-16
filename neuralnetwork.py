

"""
Biological sequences can generate feature vectors for use in machine learning progarms. Here, we predict the secondary
structure of a residue in the middle of a five-amino-acid sequence. Both the the primary and secondary structures can be
represented as code letters but then converted to zeros and ones. Then, we can pass them to the neural network.

For a five-letetr protein sequence, the vtctor will have 20 elements (one for each amino acid) for each of the five
sequence positions. The total length is 100.
"""

seqSecStrucData = [("ADTLL", "E"), ("DTLLI", "E"), ("TLLIL", "E"),
                    ("LLILG", "E"), ("LILGD", "E"), ("ILGDS", "E"),
                    ("LGDSL", "C"), ("GDSLS", "H"), ("DSLSA", "H"),
                    ("SLSAG", "H"), ("LSAGY","H"), ("SAGYR", "C"),
                    ("AGYRM", "C"), ("GYRMS", "C"), ("YRMSA", "C"), ("RMSAS","C")]

# primary structure codes
aminoAcids = "ACDEFGHIKLMNPQRSTVWY"
aaIndexDict = {}
for i, aa in enumerate(aminoAcids):
    aaIndexDict[aa] = i

# secondary structure codes
ssIndexDict = {}
ssCodes = "HCE"
for i, code in enumearte(ssCodes):
    ssIndexDict[code] = i
