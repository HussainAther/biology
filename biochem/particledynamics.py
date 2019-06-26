import numpy as np

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
    """
    Iterate through steps of the atoms in the molecule interacting
    with one another.
    """
    atoms = list(bond.keys())
    numAtoms = len(atoms)
    atomCoords = uniform(-10.0, 10.0, (numAtoms, 3))
    indices = range(numAtoms)
    n = float(numSteps)
    for step in range(numSteps):
        temp = exp(-step/n)
    for i in indices[1:]:
        atom = atoms[i]
        coords = atomCoords[i]
        velocity = np.zeros(3, float)
        for j in indices:
            if i == j:
                continue
            delta = coords - atomCoords[j]
            delta2 = delta**2
            dist2 = delta2.sum()
            bound = bondDict[atoms[j]]
            if atom in bound:
                force = bondLen - np.sqrt(dist2)
            else:
                force = 1.0 / (dist2**2)
            force = min(max(-200.0, force), 200.0)
            velocity += delta * force * temp * timeStep
            atomCoords[i] += velocity
    center = atomCoords.mean(axis=0)
    atomCoords = atomCoords - center
    return atomCoords

print(chemParticleDynamics(chemBonds))

