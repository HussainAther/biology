"""
The spectral product of A and B is the a x b two-dimensional matrix with nm
ones corresponding to all pairs of indices (ai, bj) and all remaining elements
zero. This is a psectral convolution if we take an delta-shifted peaks count
as the number of ones on the diagonal that is delta away from the main diagonal.
"""

# Create a dictionary to store each value in the spectral convolution multiset 
# along with its count.
spec = dict()
for i in A:
    for j in B:
        dist = str(i - j)
	if dist in spec:
	    spec[dist] += 1
	else:
	    spec[element] = 1

# Get the maximum multiplicity from the spectral convolution.
max_mult = max([i for i in spec.values()])

# Get the keys corresponding to the maximum multiplicities.
max_x = [item[0] for item in spec.items() if item[1] == max_mult]

