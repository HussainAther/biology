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
c_align = alignIO.read(filename, "clustal", alphabet = alpha) # from the output of a Clustalw alignment
summary_align = AlignInfo.SummaryInfo(c_align) # create the summary object from which we can determine substitutions


# Look at only amino acids with polar charged side chains by ignoring non-polar amino acids
replace_info = summary_align.replacement_dictionary(["G", "A", "V", "L", "I",
                                                    "M", "P", "F", "W", "S",
                                                    "T", "N", "Q", "Y", "C"])

