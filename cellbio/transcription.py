"""
Utilities for transcription of DNA.
"""

class Sequence:
    """
    Storing general information about sequences
    """
    TranscriptionTable = {"A":"U","T":"A","C":"G","G":"C"}
    CompTable = {"A":"T","T":"A","C":"G","G":"C"}
    def __init__(self, seqstring):
        self.seqstring=seqstring.upper()
    def restriction(self,enz):
        EnzDict={"EcoRI":"ACTGG","EcoRV":"AGTGC"}
        if EnzDict.get("EcoRI") in self.seqstring:
            return self.seqstring.count(EnzDict[enz])
        else:
            return 0
    def __getitem__(self,index):
        return self.seqstring[index]
    def __getslice__(self,low,high):
        return self.seqstring[low:high]
    def __len__(self):
        return len(self.seqstring)
    def __str__(self):
        if len(self.seqstring)>=28:
            return self.seqstring[:25]+"..."+self.seqstring[-3:]
        else:
            return self.seqstring
    def transcription(self):
        tt = ""
        for x in self.seqstring:
            if x in ’ATCG’:
                tt += self.TranscriptionTable[x]
        return tt
    def complement(self):
        tt=""
        for x in self.seqstring:
            if x in ’ATCG’:
                tt += self.CompTable[x]
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
    Reading a DNA sequence in the forward direction
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


