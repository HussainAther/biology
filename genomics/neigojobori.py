import numpy as np

"""
This method computes the numbers of synonymous and nonsynonymous substitutions and the numbers
of potentially synonymous and potentially nonsynonymous sites (Nei and Gojobori 1986).
Uses two homologous ORF sequences as input. The algorithm uses the alignment and separately
counts the number of non-synonymous and synonymous sites and differenecse between the sequences.
"""

aadict = {
"UUU":"Phe",
"UUC":"Phe",
"UUA":"Leu",
"UUG":"Leu",
"UCU":"Ser",
"UCC":"Ser",
"UCA":"Ser",
"UCG":"Ser",
"UAU":"Tyr",
"UAC":"Tyr",
"UAA":" *",
"UAG":" *",
"UGU":"Cys",
"UGC":"Cys",
"UGA":" *",
"UGG":"Trp",
"CUU":"Leu",
"CUC":"Leu",
"CUA":"Leu",
"CUG":"Leu",
"CCU":"Pro",
"CCC":"Pro",
"CCA":"Pro",
"CCG":"Pro",
"CAU":"His",
"CAC":"His",
"CAA":"Gln",
"CAG":"Gln",
"CGU":"Arg",
"CGC":"Arg",
"CGA":"Arg",
"CGG":"Arg",
"AUU":"Ile",
"AUC":"Ile",
"AUA":"Ile",
"AUG":"Met",
"ACU":"Thr",
"ACC":"Thr",
"ACA":"Thr",
"ACG":"Thr",
"AAU":"Asn",
"AAC":"Asn",
"AAA":"Lys",
"AAG":"Lys",
"AGU":"Ser",
"AGC":"Ser",
"AGA":"Arg",
"AGG":"Arg",
"GUU":"Val",
"GUC":"Val",
"GUA":"Val",
"GUG":"Val",
"GCU":"Ala",
"GCC":"Ala",
"GCA":"Ala",
"GCG":"Ala",
"GAU":"Asp",
"GAC":"Asp",
"GAA":"Glu",
"GAG":"Glu",
"GGU":"Gly",
"GGC":"Gly",
"GGA":"Gly",
"GGG":"Gly"
}

def neiGojobori(a, b):
    """
    Count non-synonymous and synonymous substitutions between DNA sequecnes a and b.
    """
    aalista = [] # list of amino acids from the a sequence
    aalistb = [] # list of amino acids from the b sequence
    aaa = "" # amino acid for a
    aab = "" # amino acid for b
    for i in range(len(a)):
        if len(aaa) == 3:
            aalista.append(aadict[aaa]) # append the amino acid to the list of amino acids in a
            aaa = ""
        else:
            
    for ind in range(len(aalista)): # for each amino acid from sequence a


    
