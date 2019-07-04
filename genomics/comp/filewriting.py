"""
For writing to files.

Draw a sequence trace using SVG polylines.
"""

# Initialize.
colors = ("black", "green", "red", "blue") 
axiscolor = "#333"
axisthickness = 2
tickcolor = "#666"
tickspacing = 100 
tickheight = 18
tickthickness = 2
maxval = 1600
leftmargin = 10
topmargin = 20
yzero = 100
hscale = 5/2 
vscale = 14

def read_data(infilname):
    """
    Read input file.
    """
    with open(infilname) as infil:
        return [eval(line) for line in infil.readlines()]
