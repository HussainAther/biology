import discrete_distribution

"""
Draw numbers from discrete probability distribution.
This could include binomial, Poisson, and negative 
binomial distributions. In biology they're used in 
gene distributions.
"""

discrete_distribution.seed()

print(sum([discrete_distribution.draw([0.2,0.2,0.3,0.3])<= for x in range(10000)])/10000.)
