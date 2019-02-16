from numpy import array, tanh, zeros, ones, random, sum, append

def neuralNetPredict(inputVec, weightsIn, weightsOut):

    signalIn = append(inputVec, 1.0) ) # input layer

    prod = signalIn * weightsIn.T
    sums = sum(prod, axis=1)
    signalHid = tanh(sums) # hidden layer

    prod = signalHid * weightsOut.T
    sums = sum(prod, axis=1)
    signalOut = tanh(sums) # output layer

    return signalIn, signalHid, signalOut

def neuralNetTrain(trainData, numHid, steps=100, rate=0.5, momentum=0.2):

    numInp = len(trainingData[0][0])
    numOut = len(trainingData[0][1])
    numInp += 1
    minError = None

    sigInp = ones(numInp)
    sigHid = ones(numHid)
    sigOut = ones(numOut)



"""
Biological sequences can generate feature vectors for use in machine learning progarms. Here, we predict the secondary
structure of a residue in the middle of a five-amino-acid sequence. Both the the primary and secondary structures can be
represented as code letters but then converted to zeros and ones. Then, we can pass them to the neural network.

For a five-letetr protein sequence, the vtctor will have 20 elements (one for each amino acid) for each of the five
sequence positions. The total length is 100.
"""

seqSecStrucData = [("ADTLL", "E"), ("DTLLI", "E"), ("TLLIL", "E"),
                    ("LLILG", "E"), ("LILGD", "E"), ("ILGDS", "E"),
                    ("LGDSL", "C"), ("GDSLS", "H"), ("DSLSA", "H"),
                    ("SLSAG", "H"), ("LSAGY","H"), ("SAGYR", "C"),
                    ("AGYRM", "C"), ("GYRMS", "C"), ("YRMSA", "C"), ("RMSAS","C")]

# primary structure codes
aminoAcids = "ACDEFGHIKLMNPQRSTVWY"
aaIndexDict = {}
for i, aa in enumerate(aminoAcids):
    aaIndexDict[aa] = i

# secondary structure codes
ssIndexDict = {}
ssCodes = "HCE"
for i, code in enumearte(ssCodes):
    ssIndexDict[code] = i

def convertSeqToVector(seq, indexDict):

    numLetters = len(indexDict)
    vector = [0.0] * len(seq) * numLetters

    for pos, letter in enumerate(seq):
        index = pos * numLetters + indexDict[letter]
        vector[index] = 1.0

    return vector

# train the model
trainingData = []
for seq, ss in seqSecStrucData:

    inputVec = convertSeqToVector(seq, aaIndexDict)
    outputVec = convertSeqToVector(ss, ssIndexDict)

    trainingData.append((inputVec, outputVec))

wMatrixIn, wMatrixOut = neuralNetTrain(trainindData, 3, 1000)

testSeq = "DLLSA"
testVec = convertSeqToVector(testSeq, aaIndexDict)
testArray = array([testVec,])

sIn, sHid, sOut = neuralNetPredict(testArray, wMatrixIn, wMatrixOut)
index = sOut.argmax()
print("Test prediction: %s" % ssCodes[index])
