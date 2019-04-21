from Bio import GenBank 
from Bio import SeqIO

"""
Utilities for parsing GenBank files.
"""

def download(accession_list):
    """
    Download and save all GenBank records in
    accession_list.
    """
    try:
        handle = GenBank.download_many(accession_list)
    except:
        print("Are you connected to the internet?")
        raise
    genbank_strings = handle.read().split("//\n") 
    for i in range(len(accession_list)):
        #Save raw file as .gb
        gb_file_name = accession_list[i] + ".gb"
        f = open(gb_file_name, "w")
        f.write(genbank_strings[i]) 
        f.write("//\n")
        f.close()

def parse(accession_list):
    """ 
    Parse all records in accession_list. 
    """
    parsed = []
    for accession_number in accession_list:
        gb_file_name = accession_number + ".gb"
        print("Parsing ... ", accession_number)
        try:
            gb_file = file(gb_file_name, "r")
        except IOError:
            print("Is the file %s downloaded?" % gb_file_name)
            raise
        gb_parsed_record = SeqIO.parse(gb_file,
                                       "genbank").next()
        gb_file.close()
        print gb_parsed_record.id 
        print gb_parsed_record.seq
        parsed.append(gb_parsed_record) 
    return parsed
