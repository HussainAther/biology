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

# Search for occurrences of a protein PROSITE patterns in the sequence

# prosite refs
sp = open(sys.argv[1])
prosite_refs = get_prosite_refs(sp)
sp.close()
# sequence
sp = open(sys.argv[1])
iterator = SProt.Iterator(sp, SProt.SequenceParser())
seq = iterator.next().seq
sp.close()

for id in prosite_refs:
    print id
    pattern=get_prosite_pattern(id)
    print pattern
    p = Pattern.compile(pattern)
    m = p.search(seq)
    print "[", m.start(), ":", m.end(), "]", seq[m.start():m.end()]
