# Program to find start codons with Position weight matrices (PWM) in the TIS
# File with a line length of 200.
# The beginning of the start codon is at position 101
# 1.11.2017

# Motives are ATG,GTG,TTG 
# ------------------------------------------------------------------------
import numpy as np
import time


# ==== constants ======================
nucDict = {'A': 0, 'C': 1, 'G': 2, 'T': 3, '\n': -1}
startCodons	= ['ATG', 'GTG', 'TTG']
sequenceFile = "TIS-Ecoli.txt"
# +++ USER CHANGEBLE +++++++++++++++++++++++++++++++++++++++++++++++++++++
pwmLength = 30 		# length of the region that is subject of the PWM
pwmStartPos = 100 	# not 101 since counting from 0 
probArr = [0.25, 0.25, 0.25, 0.25] # probability of A, C, G, T. Sum shall always be 1 (sum = 0, probability will be determined automatically)
#probArr = [0, 0, 0, 0]
seqLength = 200 	# total length of sequence including linebreake 
trainingLines_end = 722 # end number (inclusive) of the training set 
testLines_start = 0		# start number (inclusive) of the test set 
offsetToPos = 0			# number of chars that get included to the right of the PWM window
testOffset = 0			# number of chars that get included to the right of the PWM window in the test code. testOffset hast to be at least <= offsetToPos
scoreThreshold = 3.2	# score (inclusive) with witch a candidate gets counted as a valid start codon 
r = 1	# pseudo-count for PWM 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ==== Calculate probability of background model ========================================
if(sum(probArr) == 0):
	file = open(sequenceFile)
	totalCharCnt = 0
	nucleotideCnt = np.zeros((4))
	for lineIdx, line in enumerate(file):
		for charIdx, char in enumerate(line):
			if(charIdx < 100 or charIdx > 103): # Way to exclude areas from background model, for ex. the real start codons 
				totalCharCnt += 1
				nucleotideCnt[nucDict[char]] += 1
	for i in range(0, 4):
		probArr[i] = nucleotideCnt[i] / totalCharCnt
	file.close()
	print '# Background model for A,G,C,T : ', probArr


file = open(sequenceFile)
freqArr = np.zeros( (4,pwmLength + offsetToPos) ) # Nucleotide freq. 0=A, 1=C, 2=G, 3=T ; columns are positions in strings
pwmArr = np.zeros( (4,pwmLength + offsetToPos) ) # Calculated Position weight Matrix

# ==== Create frequency array, by counting nucleotide apparence at all locations ========
for lineCnt, line in enumerate(file):
	if(lineCnt <= trainingLines_end):
		for idx, char in enumerate(line):
			if(idx >= pwmStartPos - pwmLength and idx < pwmStartPos + offsetToPos): # check for the pwmLenght chars before 100 (known codon)
				if( nucDict[char] >= 0 ):
					freqArr[nucDict[char]][idx - (pwmStartPos - pwmLength + offsetToPos)] += 1 # map idx to range from 0 to n
file.close()

# ==== Calculate the PWM array with log2(freq/prob) =====================================
for rowIdx, row in enumerate(freqArr):
	tmpProb = probArr[rowIdx]
	for colIdx, cell in enumerate(row):
		pwmArr[rowIdx][colIdx] = np.log2( (cell + r) / trainingLines_end / tmpProb )

#print '=== PWM Codon Finder =============', '\n' , pwmArr, '\n', '----------------------------------'

# ==== search for possible tart codons based on pwm =====================================
file = open(sequenceFile)
codonCandidateCnt = 0
codonPwmCnt = 0
codonPwmRealMatch = 0;

for lineCnt, line in enumerate(file):
	if(lineCnt >= testLines_start):
		for i in range(0, seqLength - 2): # calculate all possible codons everywhere 
			tmpStr = line[i:i + 3]
			if( tmpStr in startCodons):
				codonCandidateCnt += 1
		# now test only the testable (30 space) and check for validity 
		for i in range(pwmLength - 1, seqLength - 3): # only go until the last 3 for possible codon 
			potantialCodon_str = line[i + 1:i + 4]
			pwm_str = line[i- (pwmLength - 1): i + 1 + testOffset]
			if( potantialCodon_str in startCodons):
				score = 0
				for charIdx, char in enumerate(pwm_str):
					score += pwmArr[nucDict[char]][charIdx]
				if(score >= scoreThreshold):
					codonPwmCnt += 1
					if(i == pwmStartPos - 1): # count matches with known start codons at pos 101(in file)
						codonPwmRealMatch += 1
file.close()

print '# Potential codons: ', codonCandidateCnt, '\n', '# CodonsIn_30PWM_area: ', codonPwmCnt, '\n', '# Match with codon: ', codonPwmRealMatch
print '# Percentage of matched codons: ', float(codonPwmRealMatch) / (sum(1 for line in open(sequenceFile)) - testLines_start), ' @', scoreThreshold, 'scoreThreshold'
print '# Percentage of false codons from the found ones: ', 1 - codonPwmRealMatch / float(codonPwmCnt)



# === On Exit ===========
exit(0)



