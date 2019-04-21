import random, sys, math

"""
"Sequence" class that represents sequences, with some fields and methods helpful
for using this class in a Gibbs sampler.
"""

class Sequence:	
	seqName = "" # name of this sequence (e.g. gene name)
	sequence = "" # nucleotide sequence
	siteScores = []	# site odds ratios = P(motif under model)/P(motif under background)
	motif = -1 # current position of motif according to Gibbs sampler
	def __init__(self, name, seq):
  		"""
		Initialize a new instance of a Sequence object, and initializes our belief of
		where the motif is to be uniform across all possible motif positions
 		"""
		self.seqName = name
		self.sequence = seq
		self.siteScores = [1 for i in range(len(seq)-motifWidth+1)]+[0 for i in range(motifWidth-1)]
		self.drawNewMotifSite()
		
	def getMotif(self, *pos):
		"""
		Return the motif of length motifWidth
		can either specify a position (e.g. getMotif(4) returns the motif at position 4)
		or if no position in specified it will return the motif at self.motif (e.g. getMotif()
		returns self.sequence[self.motif : self.motif + motifWidth]
		"""
		if pos == ():
			idx = self.motif
		else:
			idx = pos[0]
		if idx < 0 or idx > len(self.sequence)-motifWidth:
			print "Error - tried to access motif beginning at",idx,", index out of bounds."
			sys.exit(1)
		else:
			return self.sequence[idx : idx + motifWidth]
			
	def drawNewMotifSite(self):
		"""
		Randomly draws a new site for this Sequence's motif based on the current
		distribution self.siteProbs
   		"""
		tot = float(sum(self.siteScores))
		siteProbs = [x/tot for x in self.siteScores]	# normalize the siteScores
		assert abs(1.0-sum(siteProbs)) < 0.00001		# check probs sum to 1 (within margin of error)
		# draw randomly according to this distribution
		r = random.random()		# returns random uniform[0,1]
		site = 0
		cumulative = siteProbs[site]
		while cumulative < r:
			site += 1
			cumulative += siteProbs[site]
		self.motif = site		
		
	def updateSiteScores(self, wmat, background):
		"""
		Updates the odds ratios for motifs beginning at each site, where odds ratio
		= P(motif | wmat) / P(motif | background) = Pm / Pb, according to current wmat
		"""
		self.siteScores = []
		for i in range(0,len(self.sequence)):
			if i > len(self.sequence)-motifWidth:
				self.siteScores.append(0)
			else:
				Pm = 1.0		# prob under model
				Pb = 1.0		# prob under background
				for j in range(0, motifWidth):
					Pb = Pb * background[self.getMotif(i)[j]]
					Pm = Pm * wmat[j][self.getMotif(i)[j]]
				self.siteScores.append(Pm/Pb)


"""
Helper functions
"""					

def readFastaFile(filename):
   	"""
	Return a list of a Sequence objects each corresponding to 
 	a sequence in the input fasta file (filename)
	"""
	try:
		fastaLines=open(filename).readlines()
	except IOError:
		# file "filename" cannot be opened, print an error message and exit
		print "Error - could not open",filename,",please check that you entered the correct file."
		sys.exit(1)
	seqList = []
	curSeq=''
	curName=None
	for line in fastaLines:
		if '>' in line:
			if  curName !=None:
				seqList.append(Sequence(curName, curSeq))
				curName=None
				curSeq=''
			curName=line.strip()[1:]	# remove first char '>' and leading/trailing whitespace
		else:
			curSeq=curSeq+line.strip().upper()
	if curName !=None:
		seqList.append(Sequence(curName, curSeq))
	
	return seqList
	

def findSimpleBackgroundModel(sequences):
	"""
	Find a background model assuming simple 0-order model.
	Return background, a dictionary mapping the four nucleotides
	to their frequencids across all sequences.
	"""
	background = {'A': 0, 'C': 0, 'G': 0, 'T':0}
	for s in sequences:
		for nt in background:
			background[nt] = background[nt] + s.sequence.count(nt)

	# normalize
	totCounts = float(sum(background.values()))
	for nt in background:
		background[nt] = background[nt]/totCounts
	return background
	
	
def buildWeightMatrix(seqsToScore):
	"""
	Build a weight matrix from motifs in all sequences except the 
	leaveOut sequence. It includes pseudocounts at each position.
	Retuen wmat, a list of dictionary in which each position of a 
	motif is a dictionary with keys as nucleotides desrcibing the nt
	nucleotide distribution at that position.
	"""
	# initialize with pseudocounts at each position
	wmat = []
	for i in range(0, motifWidth):
		wmat.append({'A': 1, 'C': 1, 'G': 1, 'T': 1})
	# loop through all motifs, add 1 to appropriate position and nt in wmat
	for s in seqsToScore:
		for j in range(0, motifWidth):
			wmat[j][s.getMotif()[j]]+=1		
	# normalize counts
	for i in range(0, motifWidth):
		totCounts = float(sum(wmat[i].values()))
		for nt in wmat[i]:
			wmat[i][nt] = wmat[i][nt]/totCounts
				
	return wmat
	

def printWeightMatrix(wmat):
	"""
	For an input wmat weight matrix in teh format specified 
	by buildWeightMatrix(), return a human-friendly version.
	"""
	print ("Pos\tA\tC\tG\tT")
	for i in range(0,motifWidth):
		print str(i)+'\t'+str(wmat[i]['A'])+'\t'+str(wmat[i]['C'])+'\t'+str(wmat[i]['G'])+'\t'+str(wmat[i]['T'])
	
	
def calcRelEnt(wmat, background):
	"""
	Calculate the relative entropy of the weight matrix model
	wmat to the backgronud while assuming every position is 
	independent.
	"""
	relEnt = 0.0
	for pos in wmat:
		for base in pos:
			relEnt = relEnt + pos[base]*math.log(pos[base]/background[base],2)
	return relEnt
	
	
def getMotifScore(sequences, wmat, background):
	# the total score of a motif = sum of log2 (odds ratios for each sequence)
	score = 0
	for s in sequences:
		s.updateSiteScores(wmat, background)			# update with final weight matrix
		score += math.log(s.siteScores[s.motif],2)		# get score at motif
	return score
	

def printToLogo(sequences):
	# prints to the command line the motifs from each sequence in the correct format
	# for WebLogo
	for s in sequences:
		print s.getMotif()


def run(numIter):		# this is the main function

	# Get file name and motifWidth from command line
	if len(sys.argv)<=2:
		print "Error - please specify both a fasta file and motif width as inputs"
		sys.exit(1)
	else:
		fastaName=sys.argv[1]
		global motifWidth			# will make motifWidth visible to all subfunctions above
		motifWidth=int(sys.argv[2])
		
	# Read in sequences (see readFastaFile)
	sequences = readFastaFile(fastaName)
	
	# Find the background nucleotide distributions
	background = findSimpleBackgroundModel(sequences)
		
	# (STEP 1) pick starting sites at random
	# done when initializing each Sequence object in the readFastaFile() function
	
	# Repeat the following steps 2000 times
	
	for iter in range(numIter):
	
		# (STEP 2) choose sequence to leave out
		leaveOut = random.randint(0,len(sequences)-1)		# index of sequence to be left out
		seqsToScore = sequences[:]							# make list of sequences with that element
		del seqsToScore[leaveOut]							# left out
		
		# (STEP 3) make weight matrix using remaining sequences
		wmat = buildWeightMatrix(seqsToScore)
		
		# (STEP 4) update scores across all possible motif sites for left out sequence according to wmat
		sequences[leaveOut].updateSiteScores(wmat, background)
		
		# (STEP 5) draw a new site for motif in left out sequence at random according to new distribution
		sequences[leaveOut].drawNewMotifSite()
#		print calcRelEnt(wmat, background)
	
	# Print final motif matrix, its total score and its relative entropy compared to background:
#	printToLogo(sequences)
	printWeightMatrix(wmat)
	print "Motif score =",getMotifScore(sequences, wmat, background)
	print "Relative entropy =",calcRelEnt(wmat, background)
	print "Background =",background
	
	
	
run(1000)	# run main function, argument = # of iterations to run Gibbs sampler

