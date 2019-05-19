"""
Goodwin proposed a simple model for regulatory mechanisms that depend upon
feedback control. Certain metabolites may repress enzymes essential for their
own synthesis. They inhibit transcription of the molecule DNA to messenger
RNA that forms the enzyme template.
"""

def dPdt(dE, e, k, P):
    """
    According to Michaelis-Menten kinetics, degradation of product saturates  
    for large P depends upon differential of E (enzyme concentration), P
    product concentraiton, reaction rate k, and constant e. 
    """
