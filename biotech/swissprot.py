from Bio import SwissProt

# SwissProt implementation to retrieve SwissProt information

fh = open("/home/sb/bioinfo/spfile.txt")
record = SwissProt.parse(fh).next()
for att in dir(record):
    if not att.startswith(’__’):
    print(att,getattr(record,att))
