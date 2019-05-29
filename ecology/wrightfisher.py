from scipy.misc import comb

"""
Wright-Fisher model of genetic drift. (wright fisher).
The Wright-Fisher model is one of the simplest models of genetic drift in a population. 
It assumes that: (1) the generations of a population do not overlap, so that we can 
view generation change as discrete time steps (e.g., annual plants); (2) the population is 
diploid and has a finite constant size N that does not change between generations;
(3) if m alleles are present with proportion p=m2N in one generation, then each 
of the 2N chromosomes in the next generation selects the allele with probability p.

First we determine the probabiliy of a given of recessive allels in the 
first generation. We use a binomial random variable with the given parameters.
We omit the 0th term throughout the problem, as it has no contribution to the desired probability.
Then, we determine the probabiliy of a given of recessive allels in the 2nd to k-th generations.
Use the total law of probability, along with the probabilities from the previous generation.
Finally, we sum to get the desired probability.
"""

def wf(N, m, g, k):
    """
    Return the probability that in a population of N diploid individuals initially 
    possessing m copies of a dominant allele, we will observe after g generations at 
    least k copies of a recessive allele. All using the Wright-Fisher model.
    """
    p_rec = 1.0 - (m / (2.0*N))
    p = [comb(2*N, i) * ((p_rec)**i) * (1.0 - p_rec)**(2 * N-i) for i in range(1, 2 * N+1)]
    for gen in range(2,g+1):
        temp_p = []
        for j in range(1,2*N+1):
            temp_term = [comb(2*N, j)*((x/(2.0*N))**j)*(1.0-(x/(2.0*N)))**(2*N-j) for x in range(1,2*N+1)]
            temp_p.append(sum([temp_term[i]*p[i] for i in range(len(temp_term))]))
        p = temp_p
   prob = sum(p[k-1:])
   return prob
