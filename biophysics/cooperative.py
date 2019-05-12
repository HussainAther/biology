from scipy import integrate

"""
We can look at how molecules work toegether with one another to cooperate 
when forming macromolecules. We can show how the rod elasticity models of polymers
can use phenomenological parameters of the stretch (fractional change in length
of a segment), bend (how the rod becomes deformed), and twist density (torsional deformation).
We may describe worm-like chains (WLC wlc worm) for these semi-flexible polymers.
When a substrate binds to one enzymatic subunit, the rest of the subunits are 
stimulated and become active. Ligands can either have positive cooperativity, 
negative cooperativity, or non-cooperativity.
"""

def ee(T, L):
    """
    Elastic energy of a polymer using a one-parameter phenomenological model.
    It relies on temperature T for some length L. 
    """
    k = 1.3807e-23 # Boltzmann's constant (k_B) boltzmann Boltzmann
    P = .5 # Persistence (persistence) length of a polymer
    d2r2ds = lambda P*x: x**2 # d2rds2 is the end-to-end distance with some position vector r
                            # across a path s  
    return .5*k*integrate.quad(d2r2ds, 0, L)

def hill(ec10, ec90):
    """
    For a response curve coefficient, we may have input values ec10 and ec90 that produce 10% and 90% of
    the maximal response, respectively. We may then compute the Hill coefficient as a measure of ultra-
    sensitivity (how steep the response curve is).
    """
