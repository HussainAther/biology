import numpy as np

"""
Monte Carlo method randomly tests data points to approxiamte or guess a solution, thus getting a good idea
of what is going on by sampling only a few points rather than many. This approach is often coupled with trrgeting
the selection of data points towards the most promising or important results. This leads directly into the Markov chain
in which we can base a new selection on the last one. A famous example of this is the Metropolis-Hastings algorithm for
making the next guess for an effectively random walk with probabilistic selection that echoes thermodynamic energy.

Generating random numbers is key. For a simple Monte Carlo integration (determining the area bounded by some condition),
we'll use a circle.
"""

