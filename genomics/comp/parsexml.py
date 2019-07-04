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
