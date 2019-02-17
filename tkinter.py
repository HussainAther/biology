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
