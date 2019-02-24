from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO

record = SEqIO.read("https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.fna", "genbank")

"""
Create different objects directly to produce figures of a genome in linear and circular diagram
"""

