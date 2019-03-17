import numpy as np

"""
Baum-Welch algorithm is an expectation maximization algorithm used with Hidden Markov Models.
It's used in locating genes for identifying coding regions in prokaryotic DNA, analyzing
eukaroyitc seuqences up to one million base paris long, and finding copy number variations (CNVs)
among genome structure variation in humans.
"""

def bw(num_states, num_obs, observ):
    """
    Baum-Welch algorithm follows as:
    1. (expectation) forward-backward algorithm to get forward and backward probabilities
    2. create a two-dimensional array of probabilities p with n rows and k+1 columns
    3. (maximization) find the probabilities of being at a state with a time using the
    probability distribution for individual states.
    4. Repeate (iterate) expectation and maximization steps until the likelihood converges given the parameters
    """
    # Initialize the variabels we'll use
    Amat = np.ones((num_states, num_states)) / np.sum(np.ones((num_states, num_states)), 1) # transition between states matrix A
    Omat = np.ones((num_states, num_obs)) / np.sum(np.ones((num_states, num_obs)), 1) # observed elements matrix O
    theta = np.zeros((num_states, num_states, observ.size)) # likelihood
    while True: # until we converge, iterate
        oldA = Amat # for comparison, start with current matrices
        oldO = Omat
        Amat = np.ones((num_states, num_states)) # re-initialize
        Omat = np.ones((num_states, num_obs))
