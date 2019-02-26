from vpython.graph import *
from numpy import zeros

"""
Compute the population dynamics for a group of predators and prey following
the Lotka-Volterra model.
"""

Tmin = 0
Tmax = 500
y = zeros( (2), float)
Ntimes = 1000
y[0] = 2.0
y[1] = 1.3
h = (Tmax − Tmin)/Ntimes
t = Tmin

def f(t, y, F): # some function
    F[0] = 0.2∗y[0]∗(1 − (y[0]/(20.0) )) − 0.1∗y[0]∗y[1] # RHS 1st eq
    F[1] = − 0.1∗y[1] + 0.1∗y[0]∗y[1]; # RHS 2nd eq

def rk4(t, y, h, Neqs): # rk4 method
    F = zeros((Neqs), float)
    ydumb = zeros((Neqs), float)
    k1 = zeros((Neqs), float)
    k2 = zeros((Neqs), float)
    k3 = zeros((Neqs), float)
    k4 = zeros((Neqs), float)
    f(t, y, F)
    for i in range(0, Neqs):
        k1[i] = h*F[i]
        ydumb[i] = y[i] + k1[i]/2
    f(t + h/2, ydumb, F)
    for i in range(0, Neqs):
        k2[i] = h*F[i]
        ydumb[i] = y[i] + k2[i]/2
    f(t + h/2, ydumb, F)
    for i in range(0, Neqs):
        k3[i] = h*F[i]
        ydumb[i] = y[i] + k3[i]
    f(t + h, ydumb, F)
    for i in range(0, Neqs):
        k4[i] = h*F[i]
        y[i] = y[i] + (k1[i] + 2*(k2[i] + k3[i]) + k4[i]) / 6

graph1 = gdisplay(x= 0,y= 0, width = 500, height = 400,
    title = "Prey p(green) and predator P(yellow) vs time" ,
    xtitle = "t", ytitle = "P, p", xmin=0, xmax=500, ymin=0, ymax=3.5)
