import genbank 

from Bio import Seq
from Bio.Alphabet import IUPAC

"""
Genome class wraps Biopython parsing code for more effective utility.
"""

class Genome(object):
    """
    Genome - representing a genomic DNA sequence with genes
    Genome.genes[i] returns the CDS sequences for each gene i.
    """
    def __init__(self, accession_number):
        """
        Initialize by downloading from GenBank.
        """
        genbank.download([accession_number]) 
        self.parsed_genbank = genbank.parse([
                              accession_number])[0]
        self.genes = []
        self._parse_genes()
    def _parse_genes(self):
        """
        Parse out the CDS sequence for each gene.
        """
        for feature in self.parsed_genbank.features: 
            if feature.type == ’CDS’:
                #Build up a list of (start,end) tuples that will
                #be used to slice the sequence in
                #self.parsed_genbank.seq

                #Biopython locations are zero-based so can be
                #directly used in sequence splicing
                locations = []
                if len(feature.sub_features): # (4)
                    # If there are sub_features, then this gene
                    # is made up of multiple parts.  Store the
                    # start and end positins for each part.
                    for sf in feature.sub_features:
                        locations.append((sf.location.start.position,
                                        sf.location.end.position))
                else:
                    # This gene is made up of one part.  Store
                    # its start and end position.
                    locations.append((feature.location.start.position,
                                    feature.location.end.position))

                # Store the joined sequence and nucleotide
                # indices forming the CDS.
                seq = ’’ # (5)
                for begin,end in locations:
                    seq += self.parsed_genbank.seq[begin:end].tostring()
                # Reverse complement the sequence if the CDS is on
                # the minus strand
                if feature.strand == -1:  # (6)
                    seq_obj = Seq.Seq(seq,IUPAC.ambiguous_dna)
                    seq = seq_obj.reverse_complement().tostring()
                # append the gene sequence
                self.genes.append(seq) # (7)
