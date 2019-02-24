from Bio import SeqIO

"""
Identify possible genes by looking for open reading frames (ORF). We need all sixx frames for long regions
without stop codons. An ORF is a region of nucleotides with no in frame stop codons.
"""

record = SeqIO.read("https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.fna", "fasta")
table = 11 # It's a bacterial sequence so we use NCbi codon table 11 for the appropriate translation process.
min_pro_len = 100 
