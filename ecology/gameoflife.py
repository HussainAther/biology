#!/usr/bin/env python
# -*- coding: utf-8 -*- 

import vpython as vp
import numpy as np

from random import random

"""
Conway's Game of Life. A cell can either be dead (0) or alive (1).
+ If a cell is alive:
+ on the next step it will remain alive if:
+ 2 or 3 of its closer 8 neighbors are alive
+ If 3 or more neighrs are alive, the cell dies of overcrowdedness.
+ If fewer than 2 neighbors are alive, the cell dies of loneliness.
+ A dead cell will be alive if 3 of its 8 neighbors are alive.
"""

scene = vp.graph.display(width=500, height=500, title="Game of Life")
cell = np.zeros((50, 50))
cellu = np.zeros((50,50))
vp.graph.curve(pos=[(-49, -49), (-49, 49), (49, 49), (49, -49), (-49, -49)], color=color.white)
boxes = points(shape="square", size=8, color=color.cyan)

def drawcells(cell):
    """
    Draw each cell.
    """
    boxes.pos = []
    for j in range(0, 50):
        for i in range(0, 50):
            if cell[i][j] == 1:
                xx = 2*i-50
                yy = 2*j-50
                boxes.append(pos=(xx, yy))

def initial():
    """
    Initialize the board.
    """
    for j in range(20, 28):
        for i in range(20, 28):
            r = int(random()*2)
            cell[j, i] = r
    return cell

def gameoflife(cell):
    """
    Run the Game of Life simulator and check for each cell.
    """
    for i in range(1, 49):
        for j in range(1, 49):
            sum1 = cell[i−1][j−1] + cell[i][j−1] + cell[i+1][j−1] # sum 8 neighbors
            sum2 = cell[i−1][j] + cell[i+1][j] + cell[i−1][j+1] + cell[i][j+1] + cell[i+1][j+1]
            alive = sum1 + sum2
            if cell[i][j] == 1:
                if alive == 2 or alive == 3: # remain alive
                    cellu[i][j] = 1 # alive
                if alive > 3 or alive < 2: # overcrowded or loneliness
                    cellu[i][j] = 0 # rip in peace
            if cell[i][j] == 0:
                if alive == 3: # revive
                    cellu[i][j] = 1
                else:
                    cellu[i][j] = 0
    alive = 0
    return cells

temp = initial()
drawcells(temp)
while True:
    vp.rate(6)
    cell = temp
    temp = gameoflife(cell)
    drawcells(cell)
