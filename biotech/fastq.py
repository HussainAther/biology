from Bio import SeqIO

"""
Simple quality control for fastq files. Works on URLs and present directories.
"""
count = 0
for rec in SeqIO.parse("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR020/SRR020192/SRR020192.fastq", "fastq"):
    count += 1
print("%i reads" % count)
