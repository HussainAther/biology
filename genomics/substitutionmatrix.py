from Bio import AlignIO
from Bio import Alphabet
from Bio.Alphabet import IUPAC
from Bio.Align import AlignInfo

"""
Substitution matrices provide scoring terms for classifying how likely two different
residues are to substitute for each other. It's essential for sequence comparisons.
"""

filename = "http://biopython.org/DIST/docs/tutorial/examples/protein.aln"
alpha = Alphabet.Gapped(IUPAC.protein)
