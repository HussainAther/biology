import numpy as np

"""
We use a feed-backward Hidden Markov Model in predicting probabilities for how proteins fold.
"""

def forwardBackward(obs, pStart, pTrans, pEmit):
    n = len(obs)
    nStates = len(pStart)
    I = identity(nStates)
    fwd = empty([n+1, nStates])
    fwd[0] = pStart
    
    for i, val in enumerate(obs):
        fProb = np.dot(pEmit[:, val]*I, np.dot(pTrans, fwd[i]))
        fwd[i+1] = fProb / fProb.sum()

    nwd = np.ones(nStates)
    smooth = np.
