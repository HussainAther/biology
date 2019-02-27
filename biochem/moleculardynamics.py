from vpython.graph import *
import random
import numpy as np

"""
Perform a one-dimensional molecular dynamics simulation with too small a number of large steps
for just a few particles. To be realistic, we need to change the parameters and the number of random
numbers added to form the Gaussian distribution.
"""
scene = display(x=0,y=0,width=350,height=350, title="Molecular Dynamics", range=10)
sceneK = gidsplay(x=0, y=350, width=600, height=150, title="Average KE", ymin=0, ymax=5,
                    xmin=0, xmax=500, xtitle="time", ytitle"KE avg")
Kavegraph = gcurve(color=color.red)
sceneT = gidsplay(x=0, y=500, width=600, height=150, title="Average PE" ymin=-60,
                    ymin=-60, ymax=0, xmin=0, xmax=500, xtitle="time", yitle="PE avg")
Tcurve = gcurve(color=color.cyan)
Natom = 25
Nmax = 25
Tinit = 2

dens = 1
t1 = 0
x = np.zeros((Nmax), float) # x position
y = np.zeros((Nmax), float) # y position
vx = np.zeros((Nmax), float) # x velocity
vy = np.zeros((Nmax), float) # y velocity 

