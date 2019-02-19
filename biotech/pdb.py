from Bio.PDB.PDBParser import PDBParser

# PDBParser for pdf information

pdbfn = "/home/sb/bioinfo/1FAT.pdb"
parser = PDBParser(PERMISSIVE=1)
structure = parser.get_structure("1fat", pdbfn)
