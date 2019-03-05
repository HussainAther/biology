# genbank.py - utilities for downloading
# and parsing GenBank files

from Bio import GenBank #(1)
from Bio import SeqIO

def download(accession_list):
    """Download and save all GenBank records in
       accession_list.
    """
    try:
        handle = GenBank.download_many(accession_list) #(2)
    except:
        print "Are you connected to the internet?"
        raise
    genbank_strings = handle.read().split(’//\n’) #(3)
    for i in range(len(accession_list)):
        #Save raw file as .gb
        gb_file_name = accession_list[i]+’.gb’
        f = open(gb_file_name,’w’)
        f.write(genbank_strings[i]) #(4)
        f.write(’//\n’)
        f.close()

def parse(accession_list):
    """ Parse all records in accession_list. """
    parsed = []
    for accession_number in accession_list:
        gb_file_name = accession_number+’.gb’
        print ’Parsing ... ’,accession_number
        try:
            gb_file = file(gb_file_name,’r’)
        except IOError:
            print ’Is the file %s downloaded?’ % gb_file_name
            raise
        gb_parsed_record = SeqIO.parse(gb_file,
                                       "genbank").next() #(5)
        gb_file.close()
        print gb_parsed_record.id #(6)
        print gb_parsed_record.seq
        parsed.append(gb_parsed_record) #(7)
    return parsed
