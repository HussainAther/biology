import numpy as np

"""
By characterizing photon arrival frequency as a harmonic oscillator in terms
of raising and lowering (creatino and annihilation operators), we may use the 
electomagnetic field to write a Hamiltonian of a harmonic oscillator with
frequency oemga.

H = hbar * omega (a†a + 1/2)  

in which a† and a are creation and annihilation operators to connect states
with different numbers of quanta.
"""

def canda(n):
    """
    Simple calcuation of eigenstates of creation and for given numbers of 
    quanta. 
    """
    creation = np.sqrt(n+1)
    annihilation = np.sqrt(n)
    return creation, annihilation 
