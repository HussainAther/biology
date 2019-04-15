import Bio.Clustalw
import Bio.Align.AlignInfo

from Bio.Alphabet import IUPAC
from sys import *

"""
Position-specific scoring matrix
"""

if len(argv) == 2:
    threshold=40.0
else:
    threshold=argv[2]

align = Bio.Clustalw.parse_file(argv[1], alphabet=IUPAC.protein)
alig_len = align.get_alignment_length()
align_info = Bio.Align.AlignInfo.SummaryInfo(align)
ref_seq = align.get_seq_by_num(0)
pssm = align_info.pos_specific_score_matrix(ref_seq, chars_to_ignore = ["X"])
max = len(align_info.get_column(0))

print("Conservation above %d: " % threshold)
for pos in xrange(alig_len):
    for letter in pssm[pos].keys():
        percent = (pssm[pos][letter] / max) * 100.0
        if percent > threshold:
            print("%d %s %3.2f%s" % (pos, letter, percent, "%"))
