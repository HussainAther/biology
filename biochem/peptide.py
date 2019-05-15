"""
Peptides are chains of amino acids that are joined by peptide bonds. These bonds reduce the
weight of each amino acid by one H20 molecule. The result is called a residue. A Mass Spectrograph
can precisely measure the molecular weight (and charge and abundance) of any peptide chain.
Since the molecular weight of each of the possible 20 residues is known precisely, one can ask the question,
which combination of residues would give a particular weight? The problem is ambiguous for the entire molecule
"""

AminoAcid = {
    "A": "Alanine", "C": "Cysteine", "D": "Aspartic acid", "E": "Glutamic acid",
    "F": "Phenylalanine", "G": "Glycine", "H": "Histidine", "I": "Isoleucine",
    "K": "Lysine", "L": "Leucine", "M": "Methionine", "N": "Asparagine",
    "P": "Proline", "Q": "Glutamine", "R": "Arginine", "S": "Serine",
    "T": "Theronine", "V": "Valine", "W": "Tryptophan", "Y": "Tyrosine",
    "*": "STOP"
}

AminoAbbrv = {
    "A": "Ala", "C": "Cys", "D": "Asp", "E": "Glu",
    "F": "Phe", "G": "Gly", "H": "His", "I": "Ile",
    "K": "Lys", "L": "Leu", "M": "Met", "N": "Asn",
    "P": "Pro", "Q": "Gln", "R": "Arg", "S": "Ser",
    "T": "Thr", "V": "Val", "W": "Trp", "Y": "Tyr",
    "*": "STP"
}

Daltons = { 
    "A":  71, "C": 103, "D": 115, "E": 129,
    "F": 147, "G":  57, "H": 137, "I": 113,
    "K": 128, "L": 113, "M": 131, "N": 114,
    "P":  97, "Q": 128, "R": 156, "S":  87,
    "T": 101, "V":  99, "W": 186, "Y": 163
}

def predictPeptide(spectrum, prefix=""):
    """
    Predict possibly peptides from a given spectrum (list of molecular weights of a
    petide chain).
    """
    global peptideList
    if (len(prefix) == 0):
        peptideList = []
    current = sum([Daltons[res] for res in prefix])
    target = max(spectrum)
    if (current == target):
        peptideList.append(prefix)
    elif (current < target):
        for residue in Daltons.iterkeys():
            extend = prefix + residue
            suffix = [extend[i:] for i in xrange(len(extend))]
            for fragment in suffix:
                if (sum([Daltons[res] for res in fragment]) not in spectrum):
                    break
            else:
                predictPeptide(spectrum, extend)

def convolution(spectrum):
    """
    For a given noisy spectrum, we can compute the spectral convolution,
    add frequent masses above the same threshold to the spectrum,
    and use thta to infer the peptide sequence. convolve 
    """
    delta = {}
    for i in xrange(len(spectrum)-1):
        for j in xrange(i+1,len(spectrum)):
            diff = abs(spectrum[j] - spectrum[i])
            delta[diff] = delta.get(diff, 0) + 1
    return delta
