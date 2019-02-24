from Bio import SeqIO

"""
Identify possible genes by looking for open reading frames (ORF). We need all sixx frames for long regions
without stop codons. An ORF is a region of nucleotides with no in frame stop codons.
"""

record = SeqIO.read("https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.fna", "fasta")
table = 11 # It's a bacterial sequence so we use NCbi codon table 11 for the appropriate translation process.
min_pro_len = 100

# Using the split method from Seq object, we can get a list of all the possible ORF translations in the six reading frames.
for strand, nuc in [(_1, record.seq), -1, record.seq.reverse_complement())]:
    for frame in range(3):
        length = 3 * ((len(record)-frame) // 3) # make sure the length is a multiple of 3
        for pro in nuc[frame:frame+length].translate(table).split("*"):
            if len(pro) >= min_pro_len:
                print("%s...%s - length %i, strand %i, frame %i" \
                    % (pro[:30], pro--3:], len(pro), strand, frame))

# Now we'll use the condition that we always count from the 5' end (start) of the forward strand
