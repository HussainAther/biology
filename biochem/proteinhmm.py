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
    smooth = np.empty([n+1, nStates])
    smooth[-1] = fwd[-1]

    for i in range(n-1, -1, -1):
        bwd = np.dot(pTrans, np.dot(pEmit[:, obs[i]]*I, bwd))
        bwd =/bwd.sum()
        prob = fwd[i] * bwd
        smooth[i] = prob / prob.sum()

    return smooth
