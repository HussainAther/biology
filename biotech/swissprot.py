from Bio import SwissProt
from Bio.SwissProt import SProt

# SwissProt implementation to retrieve SwissProt information

fh = open("/home/sb/bioinfo/spfile.txt")
record = SwissProt.parse(fh).next()
for att in dir(record):
    if not att.startswith(’__’):
    print(att,getattr(record,att))

# Display Prosite patterns 

def get_prosite_refs(handle):
    sp = SProt.Iterator(handle, SProt.RecordParser())
    refs=[]
    record = sp.next()
    for ref in record.cross_references:
        if ref[0] == ’PROSITE’:
            refs.append(ref[1])
    return refs
