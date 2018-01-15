import numpy as np
import time
import sys

# ==== constants ================================
nucDict = {'A': 0, 'C': 1, 'G': 2, 'U': 3, '-': -1}
train_file_path = 'Daten/LSU_train.fasta'

# ==== open file and build a list of the lines =======
train_file = open(train_file_path)
train_list = train_file.readlines()
train_line_cnt = len(train_list)
train_line_lenght = len(train_list[1]) - 1 # without the \n


gap_count = np.zeros( (train_line_lenght) )			# matrix to store count of gaps per column 
nucleotide_cnt = np.zeros( (train_line_lenght, 4) )	# count of nucleotides i a column


# build the information about gaps and nucleotides at the positions 
for i in range(0, train_line_lenght):
	for j in range(1, train_line_cnt, 2):
		char = train_list[j][i]
		if( nucDict[char] >= 0 ):
			nucleotide_cnt[i][nucDict[char]] += 1
		else:
			gap_count[i] += 1


def get_emmition(start, stop):
	return 1.1

def get_transition(start, stop):
	return 1.2



# ++++ main +++++++++++++++++++++++++++++++++++++++++++++++++++
def main():
    match_list = []
    insert_list = []
    delete_list = []

    insert_list.append((0.0, 0.0)) # add initial insert state
    for i in range(0, train_line_lenght):
    	if gap_count[i]/(train_line_cnt/2) < 0.5:
    		match_list.append( (get_emmition(i, 1), 0.0 ) # 0.0 is a placeholder for the transition probability
    		insert_list.append( (0.0, 0.0) )
    	else:
    		
    		



# =========================
if __name__ == "__main__":
    main()