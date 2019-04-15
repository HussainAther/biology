import numpy as np

from math import factorial
from scipy.constants import hbar

"""
Photon statistics can be used in analyzinig how photoreceptor cells respond
to stimuli. We use the physics of random photo arrivals and the probabilistic
nature of perceptions to reflect it.

By characterizing photon arrival frequency as a harmonic oscillator in terms
of raising and lowering (creation and annihilation operators), we may use the 
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

"""
We can use the Poisson distribution to determine the probability of counting
n quanta in a certain state.
"""

def quanta(n):
    """
    Probability of counting n quanta in a given state given by the Poisson
    distribution. We write the coherent state as a superposition of states 
    with different numbers of quanta. 
    """
    c, a = canda(n) # get the creation and annihilation eigenstates
    m = np.abs(a)**2 # m is mean number of quanta
    return np.exp(-m) * m**(n) / factorial(n) 

def ham(omega, m, p, q):
    """
    Mechanistic realization of the harmonic oscillator as a mass m hangign
    from a spring with p and q momentum and position of the mass and angular
    frequency omega.
    """
    return (1/2m)*p**2 + (q**2)*(m*omega**2)/2

"""
Then we can describe the ground state as a Gaussian wavefunction.
"""

omega2 = hbar/(4*m*omega) # omega squared: zero point motion variance

q0 = (1/(2*np.pi*omega2)**(1/4)) * np.exp(-q**2/(4*omega2))  
