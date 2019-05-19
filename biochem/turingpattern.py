import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from pylab import *

"""
Alan Turing partial differential equation (PDE) models for reaction-diffusion systems 
with state variables u and v as chemical concentratinos, a b c and d parameters that 
determine the behavior of the reaction terms, with h and k constants, and Du and Dk as
diffusion constants.

dudt = a(u-h) + b(v-k) + DuΔ^2u

dudt = c(u-h) + d(v-k) + DuΔ^2u

These equations form the Turing pattern.
"""

n = 100 # grid size n x n
Dh = 1/n # spatial resolution
Dt = .02 # temporal resolution
a, b, c, d, h, k = 1, -1, 2, -1.5, 1, 1 # parameters
Du = .0001 # u diffusion constant
Dv = .0006 # v diffusion constant

# Initialize equation values
u = np.zeros([n, n])
v = np.zeros([n, n])
for x in xrange(n):
    for y in xrange(n):
        u[x, y] = 1 + uniform(-.03, .03) # small noise
        v[x, y] = 1 + uniform(-.03, .03)
nextu = np.zeros([n, n]) 
