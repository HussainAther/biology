import numpy as np
import wx
import random

"""
Grooming an ant colony.

Create a grid with randomly placed dead ants represented by blue squares and 
live ants as red circles. We create a 2-D grid with wrap-around borders like
cellular automata and compute f, the fraction of dead ants perceived by a living
ant in its neighborhood. 

Requires wxPython. 
"""

# grid specs
gridb = 640
gridh = 480 

# cell specs
hcells = 10 # horizontal cells
vcells = 10 # vertical cells
cellwidth = float(gridb)/hcells
cellheight = float(gridh/vcells
circradius = min(cellwidth, cellheight)/2 # circle radius

# ant specs
deadants = 20
liveants = 5

def paint(event):
    """
    When painting is done, set a new timer that we stop later.
    """
    panel.dc = wcPaintDC(panel)
    timer.Start(1000)

def update(event):
    """
    Bound event.
    """
    for i in range(0, vcells):
        for j in range(0, hcells):
            x = j * cellwidth
            y = i * cellheight
            if frame.world[i, j] == 1:
                col = "blue" # create blue squares
            else:
                col = "white" 
            panel.dc.setBrush(wx.Brush(col, wx.SOLID)) # solid color
            panel.dc.DrawRectangle(x, y, cellwidth, cellheight)
            xcen = x + (cellwidth/2) # x center for live ant (red circle)
            ycen = y + (cellheight/2) # y center 
            if frame.live_ant_tracker[i, j] == 2:
                col = "red" 
                panel.dc.SetBrush(wc.Brush(col, wx.SOLID))
                panel.dc.DrawCircle(xcen, ycen, circradius)
