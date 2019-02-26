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

def f( t , y, F): # some function
    F[0] = 0.2∗y[0]∗(1 − (y[0]/(20.0) )) − 0.1∗y[0]∗y[1] # RHS 1st eq
    F[1] = − 0.1∗y[1] + 0.1∗y[0]∗y[1]; # RHS 2nd eq

