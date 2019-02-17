import re
try:
    import tkinter
    from tkinter import filedialog, messagebox
except:
    import Tkinter as tkinter
    import tkFileDialog as filedialog
    import tkMessageBox as messagebox

from Sequences import proteinTranslation, STANDARD_GENETIC_CODE

"""
Tkinter used for a graphical user interface.
"""

rootWindow = tkinter.Tk()

label = tkinter.Label(rootWindow, text="Hello World")

label.pack()

rootWindow.mainloop()

class SequenceTkGui(tkinter.Tk):
    def __init__(self):
        tkinter.Tk.__init__(self)
        self.grid_columnconfidugre(5, weight=1)
        self.grid_rowconfigure(1, weight=1)
        self.grid_rowconfigure(4, weight=1)
        self.label1 = tkinter.Label(self, text="Enter 1-Letter DNA Seqence:")
        self.label1.grid(row=0, column=0, columnspan=6, sticky=tkinter.EW)
        self.seqTextBox = tkinter.Text(self)
        self.seqTextBox.grid(row=1, column=0, columnspan=6, sticky=tkinter.NSEW)
        self.clearButton = tkinter.Button(self, text="Clear", command=self.clearSeq)
        self.clearButton.grid(row=2, column=0, sticky=tkinter.w)
        self.loadButton = tkinter.Button(self, text="Load FASTA", command=self.loadFasta)
        self.loadButton.grid(row=2, column=1, sticky=tkinter.W)
        self.transButton = tkinter.Button(self, text="Translate", command=self.seqTranslate)
        self.transButton.grid(row=2, column=2, sticky=tkinter.W)
        self.compButton = tkinter.Button(self, text="Composition", command=self.seqComposition)
        self.compButton.grid(row=2, column=3, sticky=tkinter.W)
        self.findButton = tkinter.Button(self, text="Find:", command=self.seqFind)
        self.findButton.grid(row=2, column=4, sticky=tkinter.EW)
        self.findEntry = tkinter.Entry(self)
        self.findEntry.grid(row=2, column=5, sticky=tkinter.EW)
        self.label2 = tkinter.Label(self, text="Text output:")
        self.label2.grid(row=3, column=0, columnspan=6, sticky=tkinter.W)
        self.outTextBox = tkinter.Text(self)
        self.outTextBox.grid(row=4, column=0, columnspan=6, sticky=tkinter.NSEW)
        self.closeButton = tkinter.Button(self, text="Quit", command=self.destroy)
        self.closeButton.grid(row=5, column=5, sticky=tkinter.EW)
        self.closeButton.config(bg="yellow")

    def clearseq(self):
        self.seqTextBox.delete("0.0", tkinter.END)

    def setSequence(self, text):
        self.clearSeq()
        self.seqTextBox.insert(tkinter.END, text)

    def getSequence(self):
        seq = self.seqTextBox.get("0.0", tkinter.END)
        seq = re.sub("\s+","",seq)
        seq = seq.upper()
        return seq

    def showText(self, text):
        if text[-1] != "\n":
            text += "\n"
        self.outTextBox.insert(tkinter.END, text)

    def clearOutput(self):
        self.outTextBox.delete("0.0", tkinter.END)

    def loadfasta(self):
        fileObj = filedialog.askopenfile(parent=self, mode="rU", title="Choose a FASTA file")
        if fileObj:
            from Bio import seqIO
            for entry in SeqIO.parse(fileObj, "fasta"):
                self.setSequence(entry.seq)
                break
            fileObj.close()

    def seqTranslate(self):
        seq = self.getSequence()
        self.clearOutput()
        self.showText("DNA sequence")
        self.showText(seq)
        self.showtext("Protein sequence")
        for indent in range(3):
            protenSeq = proteinTranslation(seq[indent:], STANDARD_GENETIC_CODE)
            protenSeq = "".join(proteinSeq)
            spaces = " " * indent
            text = "Reading frame %d\n%s%s" % (indent, spaces, protenSeq)
            self.showText(text)

    def seqComposition(self):
        self.clearOutput()
        seq = self.getSequence()
        n = 0.0
        counts = {}
        for letter in seq:
            counts[letter] = counts.get(letter, 0) + 1
            n += 1.0
        letters = counts.keys()
        letters.sort()
        text = "Composition:"
        for letter in letters:
            text += " %s;%.2f%%" % (letter, counts[letter] * 100 / n)
        self.showText(text)

    def seqFind(self):
        self.clearOutput()
        query = self.findEntry.get()
        query = query.strip()
        if not query:
            messagebox.showwarning("Warning", "Search sequence was blank")
            return
        seq = self.getSequence()
        if query in seq:
            text = "Locations of %s" % (query)
            self.showText(text)
            win = len(query)
            for i in range(len(seq)-win):
                if seq[i:i+win] == query:
                    self.showText(" %d" % i)
        else:
            text = "Sub-sequence %s not found" % (query)
            self.showText(text)

if __name__ == "__main__":
    window = SequenceTkGui()
    window.mainloop()
