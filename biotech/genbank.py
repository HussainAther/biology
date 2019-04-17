import string

"""
Extracting the complete CDS from a GenBank entry
"""

def get_complete_cds(seqrecord):

    if string.find(seqrecord.description, "complete cds") == -1:
        return None
    for feature in seqrecord.features:
        if feature.key == "CDS":
            for qualifier in feature.qualifiers:
                if qualifier == "location":
                    location = qualifier.value
                    break
            break
    seq = seqrecord.seq
    return seq[location.start-1:location.end+1]
