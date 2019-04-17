import re

"""
Search for a specific pattern within a FASTA file.
"""

pattern = "[LIVM]{2}.RL[DE].{4}RLE" # pattern in re form
fh = open("/home/sb/bioinfo/prot.fas") # input file
fh.readline() # Discard the first line.

seq = ""
for line in fh:
    seq += line.strip()

rgx = re.compile(pattern)
result = rgx.search(seq)
patternfound = result.group()
span = result.span()
leftpos = span[0]-10
if leftpos<0:
    leftpos = 0
print(seq[leftpos:span[0]].lower()+patternfound+seq[span[1]:span[1]+10].lower())
fh.close()

