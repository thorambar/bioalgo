# Tool for local and global sequence alignment of Amino acids in fasta format 
# Thomas Pach 

import numpy as np
import time
import sys

# ==== consts and init values ================================
gl_lo_flag = True # true == local, false == global 
gap_penalty = 0.5 # linear gap penalty 
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

def get_sub_idx(str1, str2):
	try:
		retVal = sub_dict[str1, str2]
	except KeyError:
		print 'Error: non single letter code char found'
		exit(0)
	return retVal

def aligne_local():
	return 0
	# TODO

def aligne_global():
	return 0
	# TODO 

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


# build dictionary based on tuples of the two substituted acids [G, X] => -2 for ex
sub_file = open(sub_fileName)
sub_list = sub_file.readlines()
sub_list = [x.strip() for x in sub_list]
sub_dict = {}
for str in sub_list:
	sub_dict.update({(str[0], str[4]):int(str[7]+str[8])})

# open the to faster files and get the protein sequences as strings
seq_f1 = read_fasta_protein_toString(seq_f1_fileName)
seq_f2 = read_fasta_protein_toString(seq_f2_fileName)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def main():
	#gl_lo_flag, gap_penalty = init()
	exit(0)
	





# =========================
if __name__ == "__main__":
	main()