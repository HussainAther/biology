1 from Bio.Data.IUPACData
import protein_weights as protweight

# BioPython protein weight calculator

protseq = raw_input("Enter your protein sequence: ")
totalW = 0
for aa in protseq:
    totalW += protweight.get(aa.upper(),0)
totalW -= 18*(len(protseq)-1)
print("The net weight is: %s" % totalW)
