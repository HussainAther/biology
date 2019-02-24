from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO

record = SEqIO.read("https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.fna", "genbank")

"""
Create different objects directly to produce figures of a genome in linear and circular diagram
"""
# Feature set and objects
gd_feature_set = GenomeDiagram.FeatureSet()
for feature in record.features:
    if feature.type != "gene":
        # We're only looking at genes right now.
        continue
    if len(gd_feature_set) % 2 == 0:
        color = colors.blue
    else:
        color = colors.lightblue # alternate blue and light blue colors for the features
    