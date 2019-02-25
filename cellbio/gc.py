from Bio import SeqIO
from Bio.SeqUtils import GC
import pylab

"""
Calculate GC content of ls orchid fasta.
"""

url = "https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/ls_orchid.fasta"

gc_values = sorted(GC(rec.seq) for rec in SeqIO.parse(url, "fasta"))

pylab.plot(gc_values)
pylab.title("%i orchid sequences\nGC%% %0.1f to %0.1f" \
            % (len(gc_values), min(gc_values), max(gc_values)))
pylab.xlabel("Genes")
pylab.ylabel("GC%")
