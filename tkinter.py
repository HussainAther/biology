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
