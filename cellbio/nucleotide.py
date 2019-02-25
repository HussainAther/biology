from Bio import SeqIO

"""
Use dot plot to visually compare nucleotide sequences for similarity to each other.
"""

url = "https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/ls_orchid.fasta"

with open(url) as in_handle:
    record_iterator = Seq().parse(in_handle, "fasta")
    rec_one = next(record_iterator)
    rec_two = next(record_iterator)

window = 7
seq_one = str(rec_one.seq).upper()
seq_two = str(rec_two.seq).upper()
data = 
