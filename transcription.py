class Sequence:
    TranscriptionTable = {"A":"U","T":"A","C":"G","G":"C"}
    # Dictionary with the name of the restriction enzyme and
    # the recognition sequence.
    EnzDict = {"EcoRI":"GAATTC","EcoRV":"GATATC"}
    def __init__(self, seqstring):
        self.seqstring = seqstring.upper()
    def restriction(self,enz):
        try:
            return self.seqstring.count(Sequence.EnzDict[enz])
        except:
            return 0
    def transcription(self):
        tt = ""
        for x in self.seqstring:
            if x in ’ATCG’:
                tt += Sequence.TranscriptionTable[x]
        return tt

class Plasmid(Sequence):
    """
    Simple plasmid dictionary
    """
    AbResDict = {"Tet":"ctagcat","Amp":"CACTACTG"}
    def __init__(self,seqstring):
        Sequence.__init__(self,seqstring)
    def ABres(self,ab):
        if self.AbResDict[ab] in self.seqstring:
            return True
        else:
            return False

class Forward:
    """
    Reading a DNA sequence in a forward direction
    """
    def __init__(self, data):
        self.data = data
        self.index = 0
    def __iter__(self):
        return self
    def next(self):
        if self.index == len(self.data):
            raise StopIteration
        answer = self.data[self.index]
        self.index = self.index + 1
        return answer

class Reverse:
    """
    Reading a DNA sequence in reverse
    """
    def __init__(self, data):
        self.data = data
        self.index = len(data)
