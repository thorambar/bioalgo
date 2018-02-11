# Tool for local and global sequence alignment of Amino acids in fasta format 
# Thomas Pach 

import numpy as np
import time
import sys

# ==== consts and init values ================================
gl_lo_flag = True # true == local, false == global 
gap_penalty = 1 # linear gap penalty 
seq_f1_fileName = 'data/t1.fasta'
seq_f2_fileName = 'data/t2.fasta'
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

def align_local(seq_a, seq_b):
	score = 0 
	alignment = ['', '']
	len_a = len(seq_a)
	len_b = len(seq_b)
	# init dp matrix
	dp_mat = np.zeros((len_a, len_b))
	# compute score at every position 
	for i in range(1, len_a):
		for j in range(1, len_b):
			char_a = seq_a[i]
			char_b = seq_b[j]
			dp_mat[i][j] = max( (dp_mat[i-1][j-1] + get_sub_cost(char_a, char_b)), (dp_mat[i-1][j] - gap_penalty), (dp_mat[i][j-1] - gap_penalty), 0 )
	# traceback 
	return (score, [alignment[0][::-1], alignment[1][::-1]]) 
	



def align_global(seq_a, seq_b):
	score = 0;
	alignment = ['', '']
	len_a = len(seq_a)
	len_b = len(seq_b)
	# init dp matrix
	dp_mat = np.zeros((len_a+1, len_b+1))
	for i in range(0, len_a+1):
		dp_mat[i][0] = -i * gap_penalty
	for j in range(0, len_b+1):
		dp_mat[0][j] = -j * gap_penalty
	# compute score at every position 
	for i in range(1, len_a+1):
		for j in range(1, len_b+1):
			char_a = seq_a[i-1]
			char_b = seq_b[j-1]
			#print char_a, i, j
			#print i, j , char_a, char_b, max( (dp_mat[i-1][j-1] + get_sub_cost(char_a, char_b)), (dp_mat[i-1][j] - gap_penalty), (dp_mat[i][j-1] - gap_penalty) ), get_sub_cost(char_a, char_b)
			#print dp_mat
			dp_mat[i][j] = max( (dp_mat[i-1][j-1] + get_sub_cost(char_a, char_b)), (dp_mat[i-1][j] - gap_penalty), (dp_mat[i][j-1] - gap_penalty) )


	print dp_mat
	# traceback to find best alignment
	i = len_a
	j = len_b
	while i > 0 or j > 0:
		char_a = seq_a[i-1]
		char_b = seq_b[j-1]
		print char_a, i, j
		if i > 0 and j > 0 and dp_mat[i][j] == (dp_mat[i-1][j-1] + get_sub_cost(char_a, char_b)):
			# substitution 
			alignment[0] += seq_a[i-1]
			alignment[1] += seq_b[j-1]
			score += (dp_mat[i-1][j-1] + get_sub_cost(char_a, char_b))
			i -= 1
			j -= 1
		elif i > 0 and dp_mat[i][j] == (dp_mat[i-1][j] - gap_penalty):
			# gap in sequence b
			alignment[0] += seq_a[i-1]
			alignment[1] += '-'
			score += (dp_mat[i-1][j] - gap_penalty)
			i -= 1
		else:
			# gap in sequence a
			alignment[0] += '-'
			alignment[1] += seq_b[j-1]
			score += (dp_mat[i][j-1] - gap_penalty)
			j -= 1

	return (score, [alignment[0][::-1], alignment[1][::-1]]) 


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
				#str_b += '-'
				pass
		for j in range(0, step_size):
			if i+j < len_b:
				str_b += alg[1][i+j]
			else:
				#str_b += '-'
				pass
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
		print_alignment(alignment, 40)
		#print alignment[0][::-1]
	else:
		score, alignment = align_local(seq_f1, seq_f2)
		print score
		print_alignment(alignment, 60)
	exit(0)
	


# =========================
if __name__ == "__main__":
	main()