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
