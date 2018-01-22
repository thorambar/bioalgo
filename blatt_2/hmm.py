# Gene classification by use of Profile HMMs and Multiple alignment. 16.01.2018
import numpy as np
import time
import sys

# ==== constants ================================
nucDict = {'A': 0, 'C': 1, 'G': 2, 'U': 3, '-': -1}
train_file_path = 'Daten/LSU_train.fasta'
insert_threshold = 0.5
pce = [1, 1, 1, 1]
pct = [1, 1, 1]

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
	if start == stop:
		ret_arr = pce[:]
		non_gap_cnt = 0
		for i in range(0, 4):
			ret_arr[i] = nucleotide_cnt[start][i]
			non_gap_cnt += nucleotide_cnt[start][i]
		for i in range(0,4):
			ret_arr[i] = ret_arr[i] / non_gap_cnt
		return ret_arr
	else:
		ret_arr = pce[:]
		non_gap_cnt = 0
		for i in range(start, stop+1):
			for j in range(0, 4):
				ret_arr[j] = nucleotide_cnt[i][j]
				non_gap_cnt += nucleotide_cnt[i][j]
		for i in range(0, 4):
			ret_arr[i] = ret_arr[i] / non_gap_cnt
		return ret_arr




def get_transition(start, stop):
	return pct



# ++++ main +++++++++++++++++++++++++++++++++++++++++++++++++++
def main():
	match_list = []
	insert_list = []
	delete_list = []

	insert_list.append((pce, pct)) # add initial insert state
	for i in range(0, train_line_lenght):
		if gap_count[i]/(train_line_cnt/2) < insert_threshold:
			match_list.append( (get_emmition(i, i), get_transition(i, i) ) ) # pc is a placeholder for the transition probability
			insert_list.append((pce, pct)) 
			delete_list.append( (0, get_transition(i, i)) )
		else:
			#update insert
			inser_start_pos = i
			while i+1 < train_line_lenght and gap_count[i+1] > insert_threshold :
				i += 1
			insert_end_pos = i
			insert_list[-1] = (get_emmition(inser_start_pos, insert_end_pos), get_transition(i, i))


	print match_list




				  
				



# =========================
if __name__ == "__main__":
	main()