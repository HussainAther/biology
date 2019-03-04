import numpy as np

"""
This method computes the numbers of synonymous and nonsynonymous substitutions and the numbers
of potentially synonymous and potentially nonsynonymous sites (Nei and Gojobori 1986).
Uses two homologous ORF sequences as input. The algorithm uses the alignment and separately
counts the number of non-synonymous and synonymous sites and differenecse between the sequences.
"""

aadict = { # dictionary of amino acid three-codon sequence to the corresponding amino acid
"TTT":"Phe",
"TTC":"Phe",
"TTA":"Leu",
"TTG":"Leu",
"TCT":"Ser",
"TCC":"Ser",
"TCA":"Ser",
"TCG":"Ser",
"TAT":"Tyr",
"TAC":"Tyr",
"TAA":" *",
"TAG":" *",
"TGT":"Cys",
"TGC":"Cys",
"TGA":" *",
"TGG":"Trp",
"CTT":"Leu",
"CTC":"Leu",
"CTA":"Leu",
"CTG":"Leu",
"CCT":"Pro",
"CCC":"Pro",
"CCA":"Pro",
"CCG":"Pro",
"CAT":"His",
"CAC":"His",
"CAA":"Gln",
"CAG":"Gln",
"CGT":"Arg",
"CGC":"Arg",
"CGA":"Arg",
"CGG":"Arg",
"ATT":"Ile",
"ATC":"Ile",
"ATA":"Ile",
"ATG":"Met",
"ACT":"Thr",
"ACC":"Thr",
"ACA":"Thr",
"ACG":"Thr",
"AAT":"Asn",
"AAC":"Asn",
"AAA":"Lys",
"AAG":"Lys",
"AGT":"Ser",
"AGC":"Ser",
"AGA":"Arg",
"AGG":"Arg",
"GTT":"Val",
"GTC":"Val",
"GTA":"Val",
"GTG":"Val",
"GCT":"Ala",
"GCC":"Ala",
"GCA":"Ala",
"GCG":"Ala",
"GAT":"Asp",
"GAC":"Asp",
"GAA":"Glu",
"GAG":"Glu",
"GGT":"Gly",
"GGC":"Gly",
"GGA":"Gly",
"GGG":"Gly"
}

def neiGojobori(a, b):
    """
    Count non-synonymous and synonymous substitutions between DNA sequecnes a and b.
    """
    ns = 0 # non-synonomous substitutions (alters amino acid)
    ss = 0 # synonomous substitutions (doesn't alter amino acid)
    aalista = [] # list of amino acids from the a sequence
    aalistb = [] # list of amino acids from the b sequence
    aaa = "" # current amino acid for a
    aab = "" # current amino acid for b
    for i in range(len(a)):
        if len(aaa) == 3:
            aalista.append(aadict[aaa]) # append the amino acid to the list of amino acids in a
            aaa = ""
        else:
            aaa += i # add it to the current amino acid string
    for i in range(len(b)):
        if len(aab) == 3:
            aalistb.append(aadict[aab]) # append the amino acid to the list of amino acids in a
            aab = ""
        else:
            aab += i # add it to the current amino acid string
    for ind in range(len(aalista)): # for each amino acid from sequence a
        if aalista[ind] == aalistb[ind]: # the two amino acids are the exact same three strings of DNA
            pass
        else: # they aren't the exact same three strings of DNA
            if aadict[aalista[ind]] == aadict[aalistb[ind]]: # they code for the same amino acid
                ss += 1
            else: # they don't code for the same amino acid
                ns += 1
    print("Non-synonomous sbustitutions " + str(ns))
    print("Synonomous sbustitutions " + str(ss))

    print("SS p distance " + str(ss / (ss + ns)) # p distance normalizes them for comparison
    print("NS p distance " + str(ns / (ss + ns))

    

    return
