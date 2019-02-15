import discrete_distribution

# draw numbers from discrete probability distribution

discrete_distribution.seed()

print sum([discrete_distribution.draw([0.2,0.2,0.3,0.3])<= for x in range(10000)])/10000.
