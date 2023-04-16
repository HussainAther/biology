from Bio import SeqIO

INPUT_FILE = '../../samples/fasta22.fas'
OUTPUT_FILE = 'fasta22_out.fas'

def retseq(seq_fh):
    """
    Parse a fasta file and store non empty records
    into the fullseqs list.
    :seq_fh: File handle of the input seqeuce
    :return: A list with non empty sequences
    """
