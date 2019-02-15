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
    AbResDict = {"Tet":"ctagcat","Amp":"CACTACTG"}
    def __init__(self,seqstring):
        Sequence.__init__(self,seqstring)
    def ABres(self,ab):
        if self.AbResDict[ab] in self.seqstring:
            return True
        else:
            return False
