import numpy as np
import math
# ==== constants ================================
nucDict = {'A': 0, 'C': 1, 'G': 2, 'U': 3, '-': -1}
train_file_path = 'Daten/LSU_train.fasta'

# ==== open file and build a list of the lines =======
train_file = open(train_file_path)
train_list = train_file.readlines()
train_line_cnt = len(train_list) # (810)
train_line_lenght = len(train_list[1]) - 1 # without the \n (5064)

nucleotide_cnt = np.zeros( (train_line_lenght, 4) )	# count of nucleotides i a column
for i in range(0, train_line_lenght):
	for j in range(1, train_line_cnt, 2):		# only every second row contains an sequence 
		char = train_list[j][i]
		if( nucDict[char] >= 0 ):
			nucleotide_cnt[i][nucDict[char]] += 1



	tmpsum = sum(nucleotide_cnt[i])
	for j in range(0, 4):
		nucleotide_cnt[i][j] = nucleotide_cnt[i][j] / tmpsum


for i in range(0, train_line_lenght):
	entrop = 0
	for j in range(0,4):
		entrop -= nucleotide_cnt[i][j] * math.log( nucleotide_cnt[i][j]+0.00000000000000000000000000001 , 2)
	print entrop


exit(0)