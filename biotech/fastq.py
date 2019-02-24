from Bio import SeqIO

"""
Simple quality control for fastq files. Works on URLs and present directories.
"""
count = 0
ftp = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR020/SRR020192/SRR020192.fastq"

for rec in SeqIO.parse(ftp, "fastq"):
    count += 1
print("%i reads" % count)

"""
Filter for read PHREd quality 20
"""
good_reads = (rec for rec in \
              SeqIO.parse(ftp, "fastq") \
              if min(rec.letter_annotations["phred_quality"]) >= 20)
count = SeqIO.write(good_reads, "good_quality.fastq", "fastq")
print("Saved %i reads" % count)

"""
Trim off primers
"""

primer_reads = (rec for rec in \
                SeqIO.parse(ftp, "fastq") \
                if rec.seq.startswith("GATGACGGTGT"))
count = SeqIO.write(primer_reads, "with_primer.fastq", "fastq")
print("Saved %i reads" % count)
