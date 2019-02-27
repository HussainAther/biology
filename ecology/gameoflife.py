from vpython.graph import *
import random
import numpy as np

"""
Conway's Game of Life. A cell can either be dead (0) or alive (1).
+ If a cell is alive:
+ on the next step it will remain alive if:
+ 2 or 3 of its closer 8 neighbors are alive
+ If 3 or more neighrs are alive, the cell dies of overcrowdedness.
+ If fewer than 2 neighbors are alive, the cell dies of loneliness.
+ A dead cell will be alive if 3 of its 8 neighbors are alive.
"""

scene = display(width=500, height=500, title="Game of Life")
cell = np.zeros((50, 50))
cellu = zeros((50,50))
curve(pos=[(-49, -49), (-49, 49), (49, 49), (49, -49), (-49, -49)], color=color.white)
boxes = points(shape="square", size=8, color=color.cyan)

def drawcells(ce):
    boxes.pos = []
    for j in range(0, 50):
        for i in range(0, 50):
            if ce[i, j] == 1:
                xx = 2*i-50
                yy = 2*j-50
                boxes.append(pos=(xx, yy))

def initial():
    for j in range(20, 28):
        for i in range(20, 28):
            r = int(random.random()*2)
            cell[j, i] = r
    return cell

def gameoflife(cell):
    for i in range(1, 49):
        for j in range(1, 49):
            sum1 = cell[i−1,j−1]+cell[i,j−1]+cell[i+1,j−1] # sum 8 neighbors
            sum2 = cell[i−1,j]+cell[i+1,j]+cell[i−1,j+1]+cell[i,j+1]+cell[i+1,j+1]
            alive = sum1+sum2
            if cell[i,j] == 1:
                if alive == 2 or alive == 3: # remain alive
                    cellu[i,j] = 1 # alive
                if alive > 3 or alive < 2: # overcrowded or loneliness
                    cellu[i,j] = 0 # rip in peace
            if cell[i,j] == 0:
                if alive == 3: # revive
                    cellu[i,j] = 1
                else:
                    cellu[i,j]] = 0
    alive = 0
    return cells

temp = initial()
drawcells(temp)
while True:
    rate(6)
    cell = temp
    temp = gameoflife(cell)
    drawcells(cell)
