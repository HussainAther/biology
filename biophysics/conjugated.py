import numpy as np

"""
In the simplest model for a conjugated molecule, the electrons which form the pi orbitals
can sit on each carbon atom with some energy we set to zero. They hop from one atom to its
neighbors. There's one relevant electron per carbon atom. If we write the Hamiltonian for the 
electrons as a matrix, then for a ring of six carbons (benzene) we have a "hopping matrix
element" -t. It's negative as electrons lower their energy by being shared among neighboring
atoms through chemical bonds. These are tight binding models in condensed matter physics
literature and Hückel's models in the chemical literature. 
"""

# Hamiltonian for benzene
# We get this by solving Schrödinger's equation Hψ = Eψ
# to find the energy eigenstates and energy levels.
H = np.matrix([[0, -t, 0, 0, 0, -t],
	       [-t, 0, -t, 0, 0, 0],
	       [0, -t, 0, -t, 0, 0],
               [0, 0, -t, 0, -t, 0],
               [0, 0, 0, -t, 0, -t],
               [-t, 0, 0, 0, -t, 0]]) 
