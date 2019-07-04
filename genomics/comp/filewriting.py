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

def write_text(fil, x, y, txt, cls):
    """
    Write text to file.
    """
    print(" <text class='{}' x='{}' y='{}'>{}</text>".
        format(cls, x+leftmargin, y+topmargin, txt),
        file=fil)

def write_line(fil, frompos, fromval, topos, toval, color, thickness):
    """
    Write line as we want it.
    """
    print(" <line x1='{}' y1='{}' x2='{}' y2='{}'".
        format(leftmargin + int(frompos * hscale),
        int(((maxval - fromval) / vscale)) + topmargin,
        leftmargin + int(topos * hscale),
        int(((maxval - toval) / vscale)) + topmargin),
        end=' ',
        file=fil),
    print("style='stroke: {}; stroke-width: {};'/>".
        format(color, thickness),
        file=fil)

def write_axes(fil, count):
    """
    Write plot axes.
    """
    print(file=fil)
    write_line(fil, 0, 0, 0, maxval, axiscolor, axisthickness)
    write_line(fil, 0, 0, count, 0, axiscolor, axisthickness)
    print(file=fil)
    tickbottom = vscale * (-axisthickness - tickheight)
    for x in range(int(tickspacing / hscale), count, int(tickspacing // hscale)):
        write_line(fil, x, tickbottom + vscale * (1 + tickheight), x, tickbottom, tickcolor, tickthickness)
        write_text(fil, int(x * hscale), int((maxval / vscale)) + tickheight +1 5, int(x * hscale), 'tick') 

def write_bases(outfil, start, bases, positions):
    """
    Write bases.
    """
    for base, pos in zip(bases, positions): 
        write_text(outfil, int(pos * hscale - 3), int((maxval / vscale)) + 15, base, base) 
    print(file=outfil)
