import numpy as np
import wx
import random

"""
Grooming an ant colony.

Create a grid with randomly placed dead ants represented by blue squares and 
live ants as red circles. We create a 2-D grid with wrap-around borders like
cellular automata and compute f, the fraction of dead ants perceived by a living
ant in its neighborhood. 

This problem has applications in simulated annealing, scheduling problems,
vehicle routing, set problems, and nanoelectronics.

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
    Bound event. The ants should cluster dead ants toegether with the ant clustering algorithm.
    If a dead ant is found, pick it up and wander. If dead ants are found, drop the ant and wander.
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
            neighx = [x-1, x, x+1] # neighborhood for the ant  
            neighy = [y-1, y, y+1]
            dead = 0 # number of dead ants
            live = 0 # number of alive ants
            for xi in neighx:
                for yi in neighy:
                    if frame.world[xi, yi] == 1  

def randomize_worlds():
    """
    First create an empty world, initialized with 0's
    """
    frame.world = np.zeros([vcells, hcells], dtype=int)

    """
    Randomly place "dead ants" with id = 1
    """
    cell_indices = random.sample(range(hcells * vcells), deadants)
    row_indices = [idx % vcells for idx in cell_indices]
    column_indices = [idx / vcells for idx in cell_indices]
    frame.world[row_indices, column_indices] = 1

    # NOTE: Create a separate 'tracker matrix' for live ants. The reason for doing this is that if live ants step onto cells
    # containing dead ants (which they are allowed to), we want to be able to track that, as the 'world' contains information
    # only about dead ants.
    # As an alternative to maintaining a separate tracker, you can simply associate a position vector (x-y coordinates) 
    # with each ant as a property of it (in object-oriented language), and just a maintain a list of live ants.

    frame.live_ant_tracker = np.zeros([vcells, hcells], dtype=int)

    # Randomly place 'live ants' with id = 2
    cell_indices = random.sample(range(hcells * vcells), liveants)
    row_indices = [idx % vcells for idx in cell_indices]
    column_indices = [idx / vcells for idx in cell_indices]
    frame.live_ant_tracker[row_indices, column_indices] = 2

app = wx.App()
frame = wx.Frame(None, title="Act Clustering", size = (640, 500))
panel = wx.Panel(frame)
panel.SetBackgroundColour("black")

panel.Bind(wx.EVT_PAINT, paint)
timer = wx.Timer(panel)
panel.Bind(wx.EVT_TIMER, update, timer)

randomize_worlds()
        
#frame.Center()
frame.Show(True)

# A mainloop is an endless cycle that catches up all events coming up to your application.
# It is an integral part of any windows GUI application.
app.MainLoop()
