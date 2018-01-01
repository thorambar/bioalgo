# Implementation for the viterbi algorithm (casino)
import numpy as np
import time

# ==== constants ================================
file_path = 'Daten/Casino.txt'	
die_alphabet = range(1, 7)
trans_matrix = [[0.95, 0.05], [0.9, 0.1]] 		# Transition probability between loaded and fair die
e_fair = [1.0/6, 1.0/6, 1.0/6, 1.0/6, 1.0/6, 1.0/6]			# Emission probability of fair an loaded die 
e_loaded = [1.0/10, 1.0/10, 1.0/10, 1.0/10, 1.0/10, 1.0/2]
die_emissions = [e_fair, e_loaded] 		# can be selected more easily by bool value for die is loaded
die_isloaded = False
die_sequence = ''


file = open(file_path)
for line_idx, line in enumerate(file):
	if line_idx == 0:
		for char in line:
			if char != '\n':
				die_sequence += char


def viterbi(seq):
	viterbi_table = np.zeros( (3, len(seq)) )	# init viterbi table
	viterbi_table[0][0] = 1

	print viterbi_table


viterbi(die_sequence)
