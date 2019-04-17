import os
import sys

"""
Convert command line arguments to variables so that we can count kmers in them.
"""

kmer_size = int(sys.argv[1])
count_cutoff = int(sys.argv[2])
# define the function to split dna
def split_dna(dna, kmer_size):
    kmers = []
    for start in range(0,len(dna)-(kmer_size-1),1):
        kmer = dna[start:start+kmer_size]
        kmers.append(kmer)
    return kmers

# create an empty dictionary to hold the counts
kmer_counts = {}
# process each file with the right name
for file_name in os.listdir("."):
    if file_name.endswith(".dna"):
    dna_file = open(file_name)
    # process each DNA sequence in a file
    for line in dna_file:
        dna = line.rstrip("\n")
        # increase the count for each k-mer that we find
        for kmer in split_dna(dna, kmer_size):
            current_count = kmer_counts.get(kmer, 0)
            new_count = current_count + 1
            kmer_counts[kmer] = new_count

# print k-mers whose counts are above the cutoff
for kmer, count in kmer_counts.items():
    if count > count_cutoff:
    print(kmer + " : " + str(count))
