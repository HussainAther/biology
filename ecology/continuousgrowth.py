import numpy as np

"""
Among continuous growth models of population, we can look at 
single-species models that are relevant to lab studies. But in the real world
we can reflect multiple types of phenomena and events that occur.

We'll define the rate of change in this simple context

dNdt = births - deaths + migration

and get a conservation equation (if we choose to ignore migration and assume
the births and deaths are proportional to N)

dNdt = bN - dN 

and solve for this

N(t) = N0 * exp((b-d)t)
"""

def N(t):
    """
    Continuous growth model as described above under some conditions of birth  
    and death.
    """
    b = 5
    d =  4 
    N0 = np.exp((b-d)*t)
    return N0 * np.exp((b-d)*t)
