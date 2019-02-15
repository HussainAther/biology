class Sequence:
    TranscriptionTable = {"A":"U","T":"A","C":"G","G":"C"}
    def __init__(self, seqstring):
        self.seqstring=seqstring.upper()
    def transcription(self):
        tt = ""
        for x in self.seqstring:
            if x in ’ATCG’:
                tt += self.TranscriptionTable[x]
        return tt
