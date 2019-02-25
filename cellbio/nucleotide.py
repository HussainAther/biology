from Bio import SeqIO
import pylab

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
data = [[(seq_one[i:i + window] != seq_two[j:j + windpw])
        for j in range(len(seq_one) - window)]
        for i in range(len(seq_two) - window)]

pylab.gray()
pylab.imshow(data)
pylab.xlabel("%s (length %i bp)" % (rec_one.id, len(rec_one)))
pylab.ylabel("%s (length %i bp)" % (rec_two.id, len(rec_two)))
pylabe.title("Dot plot using window size %i\n(allowing no mis-matches)" % window)
pylab.show()
