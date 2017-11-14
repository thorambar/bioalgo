# Program to find start codons with Position weight matrices (PWM) in the TIS
# File with a line length of 200.
# The beginning of the start codon is at position 101
# 1.11.2017

# Motives are ATG,GTG,TTG 
# ------------------------------------------------------------------------
import numpy as np


# ==== constants ======================
startCodons =['ATG', 'GTG', 'TTG']
sequenceFile = "TIS-Ecoli.txt"
pwmLength = 30 
pwmStartPos = 100 # not 101 since counting from 0 
prob = 0.25
seqLength = 200
exampleLines = 722
r = 1
nucDict = {'A': 0, 'C': 1, 'G': 2, 'T': 3, '\n': -1}


file = open(sequenceFile)
freqArr = np.zeros( (4,pwmLength) ) # Nucleotide freq. 0=A, 1=C, 2=G, 3=T ; columns are positions in strings
pwmArr = np.zeros( (4,pwmLength) ) # Calculated Position weight Matrix

# Create frequency array, by counting nucleotide apparence at all locations
for line in file:
	for idx, char in enumerate(line):
		if(idx >= pwmStartPos - pwmLength and idx < pwmStartPos): # check for the pwmLenght chars before 100 (known codon)
			if( nucDict[char] >= 0 ):
				freqArr[nucDict[char]][idx - (pwmStartPos - pwmLength)] += 1 # map idx to range from 0 to n

# Calculate the PWM array with log2(freq/prob)
for rowIdx, row in enumerate(freqArr):
	for colIdx, cell in enumerate(row):
		pwmArr[rowIdx][colIdx] = np.log2( (cell+1) / exampleLines / prob ) # TODO: Divide by 0 encountered vielleicht nicht schlimm 

print pwmArr






# === On Exit ===========
exit(0)



