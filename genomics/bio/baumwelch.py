import numpy as np

"""
Baum-Welch algorithm is an expectation maximization algorithm used with Hidden Markov Models.
It's used in locating genes for identifying coding regions in prokaryotic DNA, analyzing
eukaroyitc seuqences up to one million base paris long, and finding copy number variations (CNVs)
among genome structure variation in humans.

Going to implement this using the feedback algorithm of each individual step following equations for
determining the various values in python code.
"""

def fb(Amat, Omat, observ):
    """
    Feedback algorithm among our matrices of transition and observation.
    """
    k = observ.size # size of waht we observe
    (n,m) = Omat.shape # shape for future reference
    probmat = np.zeros((n,k)) # probability matrix for the transitions using the appropriate shape
    fw = np.zeros((n,k+1)) # forward step
    bw = np.zeros((n,k+1)) # backward step
    fw[:, 0] = 1.0/n # forward step initialized
    for obsind in xrange(k): # for the observation indices
        frowvec = np.matrix(fw[:,obsind]) # convert to matrix and save as a vector the current row
        fw[:, obsind+1] = frowvec * np.matrix(Amat) * np.matrix(np.diag(Omat[:,observ[obsind]])) # find the probability for each state and time
        fw[:, obsind+1] = fw[:,obsind+1]/np.sum(fw[:,obsind+1]) # normalize with the sum
    bw[:,-1] = 1.0 # backward step initialized. same as with forward.
    for obsind in xrange(k, 0, -1): # observation index
        bcolvec = np.matrix(bw[:,obsind]).transpose()
        bw[:, obsind-1] = (np.matrix(Amat) * np.matrix(np.diag(Omat[:,observ[obsind-1]])) * bcolvec).transpose()
        bw[:,obsind-1] = bw[:,obsind-1]/np.sum(bw[:,obsind-1])
    probmat = np.array(fw)*np.array(bw)
    probmat = probmat/np.sum(probmat, 0)
    return probmat, fw, bw

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
        p, fo, bw = fb( oldA, oldO, observ) # expectation
        for a in xrange(num_states):
            for b in xrange(num_states):
                for t in xrange(observ.size):
                    theta[a, b, t] = fo[a,t] * bw[b, t+1] * oldA[a, b] * oldO[b, observ[t]]
        for a in xrange(num_states):
            for b in xrange(num_states):
                for t in xrange(observ.size):
                    

