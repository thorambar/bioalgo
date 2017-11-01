# Program to find start codons with Position weight matrices (PWM) in the TIS
# File with a line length of 200.
# The beginning of the start codon is at position 101
# By Thomas Pach 1.11.2017

# ==== constants ======================
sequenceFile = "TIS-Ecoli.txt"
pwmLength = 30 
r = 1



file = open(sequenceFile)
lineCount = 0
for line in file:
	lineCount += 1
	
print "# File stats: ", lineCount




