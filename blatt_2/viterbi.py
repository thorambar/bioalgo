# Implementation for the viterbi algorithm (casino)
import numpy as np
import time
import sys
# ==== constants ================================
file_path = 'Daten/Casino.txt'	
trans_matrix = [[0.95, 0.05], [0.9, 0.1]] 		# Transition probability between loaded and fair die
e_fair = [1.0/6, 1.0/6, 1.0/6, 1.0/6, 1.0/6, 1.0/6]			# Emission probability of fair an loaded die 
e_loaded = [1.0/10, 1.0/10, 1.0/10, 1.0/10, 1.0/10, 1.0/2]
die_emissions = [e_fair, e_loaded] 		# can be selected more easily by bool value for die is loaded
die_sequence = ''


file = open(file_path)
for line_idx, line in enumerate(file):
	if line_idx == 0:
		for char in line:
			if char != '\n':
				die_sequence += char


def viterbi(seq):
	viterbi_table = np.zeros( (2, len(seq)+1) )	# init viterbi table
	viterbi_table[0][0] = 1.0
	viterbi_table[1][0] = 1.0

	for idx, char in enumerate(die_sequence):
		if idx < len(seq):
			viterbi_table[0][idx+1] = die_emissions[0][int(char)-1] * max( (viterbi_table[0][idx] * trans_matrix[0][0]), (viterbi_table[1][idx] * trans_matrix[1][1]) )
			viterbi_table[1][idx+1] = die_emissions[1][int(char)-1] * max( (viterbi_table[0][idx] * trans_matrix[0][1]), (viterbi_table[1][idx] * trans_matrix[1][0]) )


	if viterbi_table[0][len(seq)] > viterbi_table[1][len(seq)]:
		end_state = 0
	else:
		end_state = 1
	state = end_state
	output_stack = []
	for idx in range(len(seq)+1, 1, -1):
		if state == 0:
			if viterbi_table[0][idx-1] * trans_matrix[0][0] > viterbi_table[1][idx-1] * trans_matrix[1][1]:
				state = 0
				output_stack.append('F')
			else:
				state = 1
				output_stack.append('L')
		else:
			if viterbi_table[1][idx-1] * trans_matrix[1][0] > viterbi_table[0][idx-1] * trans_matrix[0][1]:
				state = 1
				output_stack.append('L')
			else:
				state = 0
				output_stack.append('F')

	output_stack.reverse()
	for idx, i in enumerate(output_stack):
		sys.stdout.write(i)



	


viterbi(die_sequence)
exit(0)