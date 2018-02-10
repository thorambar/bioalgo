# Tool for local and global sequence alignment of Amino acids in fasta format 
# Thomas Pach 

import numpy as np
import time
import sys

# ==== consts and init values ================================
gl_lo_flag = True # true == local, false == global 
gap_penalty = 0.5 # linear gap penalty 
seq_f1_fileName = 'data/TPA_HUMAN'
seq_f2_fileName = 'data/TPA_PIG'
sub_fileName = 'data/BLOSUM62'


# build dictionary based on tuples of the two substituted acids [G, X] => -2 for ex
sub_file = open(sub_fileName)
sub_list = sub_file.readlines()
sub_list = [x.strip() for x in sub_list]
sub_dict = {}
for str in sub_list:
	sub_dict.update({(str[0], str[4]):int(str[7]+str[8])})
print sub_dict[('S', 'W')]





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



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def main():
	#gl_lo_flag, gap_penalty = init()
	print ''





# =========================
if __name__ == "__main__":
	main()