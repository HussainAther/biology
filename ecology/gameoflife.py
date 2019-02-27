from vpython.graph import *
import random 

"""
Conway's Game of Life. A cell can either be dead (0) or alive (1).
+ If a cell is alive:
+ on the next step it will remain alive if:
+ 2 or 3 of its closer 8 neighbors are alive
+ If 3 or more neighrs are alive, the cell dies of overcrowdedness.
+ If fewer than 2 neighbors are alive, the cell dies of loneliness.
+ A dead cell will be alive if 3 of its 8 neighbors are alive.
"""
