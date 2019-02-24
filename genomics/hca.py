from Bio import Cluster

"""
Hiercarchical clustering algorithm using BioPython. We use single linkage clustering for genes
and maximum linkage clustering for experimental conditions. As the Euclidean distance is
"""

with open("https://raw.githubusercontent.com/biopython/biopython/master/Tests/Cluster/cyano.txt" as handle):
    record = Cluster.read(handle)

genetree = record.treecluster(method = "s") # s means pairwise single-linkage clustering
