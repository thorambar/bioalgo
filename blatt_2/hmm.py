# Gene classification by use of Profile HMMs and Multiple alignment. 16.01.2018
import numpy as np
import time
import sys

# ==== constants ================================
nucDict = {'A': 0, 'C': 1, 'G': 2, 'U': 3, '-': -1}
MATCH = 0
INSERT = 1
DELETE = 2
train_file_path = 'Daten/LSU_train.fasta'
insert_threshold = 0.5 # when there are more than x percent gaps, its an insert  
pce = [1, 1, 1, 1] # A, C, G, T
pct = [1, 1, 1]	# Match, Insert, Delete 

# ==== open file and build a list of the lines =======
train_file = open(train_file_path)
train_list = train_file.readlines()
train_line_cnt = len(train_list) # (810)
train_line_lenght = len(train_list[1]) - 1 # without the \n (5064)

gap_count = np.zeros( (train_line_lenght) )			# matrix to store count of gaps per column 
nucleotide_cnt = np.zeros( (train_line_lenght, 4) )	# count of nucleotides i a column


# build the information about gaps and nucleotides at the positions 
for i in range(0, train_line_lenght):
	for j in range(1, train_line_cnt, 2):		# only every second row contains an sequence 
		char = train_list[j][i]
		if( nucDict[char] >= 0 ):
			nucleotide_cnt[i][nucDict[char]] += 1
		else:
			gap_count[i] += 1


def gwc(pos):
	# get weighted gap count 
	if gap_count[pos] / (train_line_cnt/2) < insert_threshold:
		return False
	else:
		return True

def get_emission(start, stop):
	# computes the emission probability by counting all base occurrences in a column 
	# for columns from start to stop point 
	ret_arr = pce[:]
	non_gap_cnt = 0
	for i in range(start, stop+1):
		for j in range(0, 4):
			ret_arr[j] += nucleotide_cnt[i][j]
			non_gap_cnt += nucleotide_cnt[i][j]
	for i in range(0, 4):
		ret_arr[i] = ret_arr[i] / non_gap_cnt
	return ret_arr




def get_transition(start, state):
	ret_arr = pct[:]
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	if state == MATCH:
		for i in range(1, train_line_cnt, 2): 	# iterate over every other row 
			while start+1 < train_line_lenght:
				if start+1 < train_line_cnt and train_list[i][start+1] > 0 and not gwc(start+1):
					# its a match
					ret_arr[0] += 1
					break	
				elif start+1 < train_line_cnt and not gwc(start+1) and train_list[i][start+1] < 0:
					# its a delete
					ret_arr[2] += 1
					break
				elif start+1 < train_line_cnt and gwc(start+1) and train_list[i][start+1] > 0:
					# its a insert	
					ret_arr[1] += 1
					break
				start += 1
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	elif state == INSERT:
		return ret_arr
			
	return ret_arr



# ++++ main +++++++++++++++++++++++++++++++++++++++++++++++++++
def main():
	match_list = []
	insert_list = []
	delete_list = []

	insert_list.append((pce, pct)) 												# add initial insert state
	for i in range(0, train_line_lenght): 										# iterate over chars in one line
		if gap_count[i]/(train_line_cnt/2) < insert_threshold:
			match_list.append( (get_emission(i, i), get_transition(i, MATCH) ) ) 
			insert_list.append((pce, pct)) 
			delete_list.append( (0, get_transition(i, DELETE)) )
		else:																	#update insert
			inser_start_pos = i
			while i+1 < train_line_lenght and gap_count[i+1] > insert_threshold :
				i += 1
			insert_end_pos = i
			insert_list[-1] = (get_emission(inser_start_pos, insert_end_pos), get_transition(i, INSERT))



	




				  
				



# =========================
if __name__ == "__main__":
	main()