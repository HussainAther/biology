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


