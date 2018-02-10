# Tool for local and global sequence alignment of Amino acids in fasta format 
# Thomas Pach 

import numpy as np
import time
import sys

# ==== consts and init values ================================
gl_lo_flag = True # true == local, false == global 
gap_penalty = 1 # linear gap penalty 
seq_f1_fileName = 'data/TPA_HUMAN.fasta'
seq_f2_fileName = 'data/TPA_PIG.fasta'
sub_fileName = 'data/BLOSUM62'

def read_fasta_protein_toString(path):
	file = open(path)
	seq = ''
	for line in file:
		line = line.strip()
		if line.startswith('>'):
			pass
		else:
			seq += line
	return seq

def get_sub_cost(str1, str2):
	try:
		retVal = sub_dict[str1, str2]
	except KeyError:
		print 'Error: non single letter code char found'
		exit(0)
	return retVal

def align_local():
	return 0
	# TODO

def align_global(seq_a, seq_b):
	score = 0;
	alignment = ['', '']
	len_a = len(seq_a)
	len_b = len(seq_b)
	# init dp matrix
	dp_mat = np.zeros((len_a, len_b))
	for i in range(0, len_a):
		dp_mat[i][0] = -i * gap_penalty
	for j in range(0, len_b):
		dp_mat[0][j] = -j * gap_penalty
	# compute score at every position 
	for i in range(1, len_a):
		for j in range(1, len_b):
			char_a = seq_a[i]
			char_b = seq_b[j]
			dp_mat[i][j] = max( (dp_mat[i-1][j-1] + get_sub_cost(char_a, char_b)), (dp_mat[i-1][j] - gap_penalty), (dp_mat[i][j-1] - gap_penalty) )
	# traceback to find best alignment
	for i in range(len_a-1, -1, -1):
		opt_j = 0
		opt = 0
		opt_idx = 0
		for j in range(len_a-1, -1, -1):
			tmp_cost_list = [dp_mat[i-1][j-1], dp_mat[i-1][j], dp_mat[i][j-1]]
			tmp_max = max(tmp_cost_list)
			tmp_idx = tmp_cost_list.index(tmp_max)
			if tmp_max > opt:
				opt = tmp_max
				opt_j = j
				opt_idx = tmp_idx
		score += opt
		if opt_idx == 0:
			# substitution 
			alignment[0] += seq_a[i]
			alignment[1] += seq_b[opt_j]
		elif opt_idx == 1:
			# gap in sequence a
			alignment[0] += '-'
			alignment[1] += seq_b[opt_j]
		else:
			# gap in sequence b
			alignment[0] += seq_a[i]
			alignment[1] += '-'
	return (score, alignment) 


def init():
	tmp = raw_input('Press [1] for local alignment or [2] for global alignment: ')
	try: 
		int_flag = int(tmp)
	except ValueError:
		print 'Wrong input'
		exit(0)
	if(int_flag ==  2):
		flag = False
	
	tmp = raw_input('pleas enter gap penalty: ')
	try: 
		print float(tmp)
		gpen = float(tmp)
	except ValueError:
		print 'Wrong input'
		exit(0)

	return (flag, gpen)

def print_alignment(alg, step_size):
	len_a = len(alg[0])
	len_b = len(alg[1])
	i = 0
	while i < max(len_a, len_b):
		str_a = ''
		str_b = ''
		for j in range(0, step_size):
			if i+j < len_a:
				str_a += alg[0][i+j]
			else:
				str_b += '-'
		for j in range(0, step_size):
			if i+j < len_a:
				str_b += alg[1][i+j]
			else:
				str_b += '-'
		print str_a
		print str_b
		print ' '
		print ' '
		i += step_size





# build dictionary based on tuples of the two substituted acids [G, X] => -2 for ex
sub_file = open(sub_fileName)
sub_list = sub_file.readlines()
sub_list = [x.strip() for x in sub_list]
sub_dict = {}
for str in sub_list:
	sub_dict.update({(str[0], str[4]):int(str[7]+str[8])})



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def main():
	# open the to faster files and get the protein sequences as strings
	seq_f1 = read_fasta_protein_toString(seq_f1_fileName)
	seq_f2 = read_fasta_protein_toString(seq_f2_fileName)
	
	#gl_lo_flag, gap_penalty = init()
	gl_lo_flag = False

	if gl_lo_flag == False:
		score, alignment = align_global(seq_f1, seq_f2)
		print score
		print_alignment(alignment, 60)
		#print alignment[1]
	else:
		score, alignment = align_local(seq_f1, seq_f2)
		print score
		print alignment
	exit(0)
	


# =========================
if __name__ == "__main__":
	main()