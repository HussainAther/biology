import numpy as np

"""
We use a feed-backward Hidden Markov Model (HMM) in predicting probabilities for how proteins fold.
"""

def forwardBackward(obs, pStart, pTrans, pEmit):
    """
    Move forward or backward in the Hidden Markov Model.
    """
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

expTypes = ["-", "*"]
aaTypes = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

nExp = len(expTypes)
nAmino = len(aaTypes)
nStates = nExp * nAmino

# initialize probabilities
pStart = np.zeros(nStates, float)
pTrans = np.zeros((nStates, nStates), float)
pEmit = np.zeros((nSt, nAmino), float)

indexDict = {}
stateDict= {}
index = 0

# calculate and store probabilities
for exposure in expTypes:
    for aminoAcid in aaTypes:
        stateKey = (exposure, aminoAcid)
        indexDict[stateKey] = index
        stateDict[index] = stateKey
        index += 1
