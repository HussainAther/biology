from numpy import zeros, sqrt

"""
In integrating the various ways particles interact with one another, we can determine solutions that
hold true for various folding of proteins and create models of other biochemical structures. We can model a simple
bond between two atoms with a repulsive force. The test molecule is acetamoniophen, a common pharmaceutical paracetamol.
"""

chemBonds = {"H1": ["01",], "01" : ["H1", "C1"], "C1" : ["01", "C2", "C6"],
            "C2": ["C1", "H2", "C3"], "H2" : ["C2"], "C3" ["C2", "H3", "C4"],
            "H3" : ["C3"], "C4" : ["C3", "N7", "C5"], "C5" : ["C4", "H5", "C6"],
            "H5" : ["C5"], "C6" : ["C5", "H6", "C1"], "H6" : ["C6"],
            "N7" : ["C4", "H7", "C8"], "H7" ["N7"], "C8" : ["N7", "08", "C9"],
            "08" : ["C8"], "C9" : ["C8", "H9a", "H9b", "H9c"], "H9a" : ["C9"],
            "H9b" : ["C9"], "H9c" : ["C9"]
}

def chemParticleDynamics(bondDict, numSteps=5000, bondLen=1.0, timeStep=0.01):
    atoms = list(bond.keys())
    numAtoms = len(atoms)
    atomCoords = uniform(-10.0, 10.0, (numAtoms, 3))
