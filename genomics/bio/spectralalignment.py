

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

