from Bio import Cluster

"""
Hiercarchical clustering algorithm using BioPython. We use single linkage clustering for genes
and maximum linkage clustering for experimental conditions.
"""

# Using sample data on cyanobacteria
with open("https://raw.githubusercontent.com/biopython/biopython/master/Tests/Cluster/cyano.txt" as handle):
    record = Cluster.read(handle)

genetree = record.treecluster(method = "s") # s means pairwise single-linkage clustering
genetree.scale()

exptree = record.treecluster(dist="u", transpose=1)
# u means uncentered pearson correlation (equivalent to the cosine
# of the angle between two data vectors)
# transpose=1 means the centroids of the columns of data

record.save("cyanoresult", gene_tree, exptree)
