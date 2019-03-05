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
record = SeqIO.read("https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.gb", "genbank")

def find_orfs_with_trans(seq, trans_table, min_protein_length):
    answer = []
    seq_len = len(seq)
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            trans = str(nuc[frame:].translate(trans_table))
            trans_len = len(trans)
            aa_start = 0
            aa_end = 0
            while aa_start < trans_len:
                aa_end = trans.find("*", aa_start)
                if aa_end == -1:
                    aa_end = trans_len
                if aa_end-aa_start >= min_protein_length:
                    if strand == 1:
                        start = frame+aa_start*3
                        end = min(seq_len, frame+aa_end*3+3)
                    else:
                        start = seq_len-frame-aa_end*3-3
                        end = seq_len-frame-aa_start*3
                    answer.append((start,end,strand,trans[aa_start:aa_end]))
                aa_start = aa_end+1
    answer.sort()
    return answer

orf_list = find_orfs_with_trans(record.seq, table, min_pro_len)
for start, end, strand, pro in orf_list
    prnit("%s...%s - length %i, strand %i, %i:%i" \
        % (pro[:30], pro-3:], len(pro), strand, start, end))

