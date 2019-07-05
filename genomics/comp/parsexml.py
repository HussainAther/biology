import xml.parsers.expat

"""
Get content of several XML tags.
"""

class GenomeBinomialParser:
    ChunkSize = 4000 # number of bytes to read at a time
    def __init__(self):
        self.parser = xml.parsers.expat.ParserCreate()
        self.parser.pbuffer_text = True
        self.parser.pStartElementHandler = self.start_element
        self.parser.pEndElementHandler = self.end_element
        self.parser.pCharacterDataHandler = self.char_text

     def start_element(self, name, attrs):
         """
         Record name of start tag just encountered.
         """
         self.current_tag = name # ignore attrs
    
     def end_element(self, name):
         """
         Stop parsing when species name has been encountered.
         """
         if name == 'BinomialOrgName_species':
             raise StopIteration
    
     def char_text(self, text):
         """
         If current tag is for genus or species, record its content.
         """
         if not text.isspace():
             # checking for whitespace because this function gets called
             # for both tag content and whitespace between the end of
             # one tag and the beginning of the next
             if self.current_tag == 'BinomialOrgName_genus':
                 self.genus = text.strip()
             elif self.current_tag == 'BinomialOrgName_species':
                 self.species = text.strip()

    def find_binomial(self, filename):
        """
        Return (genus, species) as found in the XML genome file named filename.
        """
        # read only ChunkSize characters at a time in case we find it relatively early
        # and to not take up an enormous amount of memory with file contents
        with open(filename) as file:
            try:
                while True:
                    self.parser.Parse(file.read(self.ChunkSize))
            except StopIteration:
                pass
        return self.genus, self.species
