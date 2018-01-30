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
#train_file_path = 'Daten/my_test.fasta'
test_file_path = 'Daten/LSU_short_test.fasta'
insert_threshold = 0.5 # when there are more than x percent gaps, its an insert  
pce = [1, 1, 1, 1] # A, C, G, T
pct = [1, 1, 1]	# Match, Insert, Delete 
match_list = []
insert_list = []
delete_list = []

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
	# get weighted gap count (returns true if the column is a gap)
	if gap_count[pos] / (float(train_line_cnt)/2) < insert_threshold:
		return False
	else:
		return True


def get_emission(start, stop):
	# computes the emission probability by counting all base occurrences in a column 
	# for columns from start to stop point 
	ret_arr = pce[:]
	non_gap_cnt = sum(pce)
	for i in range(start, stop+1):
		for j in range(0, 4):
			ret_arr[j] += nucleotide_cnt[i][j]
			non_gap_cnt += nucleotide_cnt[i][j]
	for i in range(0, 4):
		ret_arr[i] = ret_arr[i] / non_gap_cnt
	return ret_arr



def get_transition(start, state):
	ret_arr = pct[:]

	#print 'trans ----------'

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	if state == MATCH:
		#print 'Match'
		for i in range(1, train_line_cnt, 2): 	# iterate over every other row 
			#print '------------------------- i: ', i
			while start+1 < train_line_lenght:
				#print 'start: ', start
				if start+1 < train_line_cnt and train_list[i][start+1] > 0 and not gwc(start+1):
					# its a match
					ret_arr[0] += 1
					#print 'Match'
					break	
				elif start+1 < train_line_cnt and not gwc(start+1) and train_list[i][start+1] < 0:
					# its a delete
					ret_arr[2] += 1
					#print 'delete'
					break
				elif start+1 < train_line_cnt and gwc(start+1) and train_list[i][start+1] > 0:
					# its a insert	
					ret_arr[1] += 1
					#print 'insert'
					break
				start += 1
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	elif state == INSERT:
		#print 'Insert'
		for i in range(1, train_line_cnt, 2):
			while start+1 < train_line_lenght:
				if not gwc(start+1): # when its a non gap state 
					if train_list[i][start+1] > 0:
						# its a match
						ret_arr[0] += 1
						#print 'Match', start
						break
					else:
						# its a delete 
						ret_arr[2] += 1
						#print 'delete', start
						break
				elif gwc(start+1) and train_list[i][start+1] > 0:
					# its a insert 
					#print 'insert', start
					ret_arr[1] += 1 # no break, since inserts of multiple nucleotides are more common
				start += 1
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	elif state == DELETE:
		#print 'Delete'
		for i in range(1, train_line_cnt, 2):
			while start+1 < train_line_lenght:
				if not gwc(start +1):
					if train_list[i][start+1] < 0:
						# its an delete
						ret_arr[2] += 1
						#print 'delete', start
						break
					if train_list[i][start+1] > 0:
						# its a match
						ret_arr[0] += 1
						#print 'Match', start
						break
				elif gwc(start+1) and train_list[i][start+1] > 0:
					# its an insert 
					ret_arr[1] += 1
					#print 'Match', start
					break
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	tmp_sum = sum(ret_arr) 
	#print ret_arr
	for i in range(0, 3):
		ret_arr[i] = float(ret_arr[i]) / tmp_sum # normalize probability
	#print ret_arr, sum(ret_arr)
	return ret_arr


def viterbi(string):
	current_state_idx = [0, 0, 0] # list containing the index of [Match, Insert, Delete] node in corresponding lists
	current_state = 0 # type of state (0,1,2) in which we currently are
	current_state_trans = (0,0,0)
	path_string = ''
	path_likelihood = 0.0
	
	# init with highest likelihood (only match, insert possible?) 
	if ( match_list[0][0][nucDict[string[0]]] > insert_list[0][0][nucDict[string[0]]] ):
		# if start is match
		path_likelihood = match_list[0][0][nucDict[string[0]]]
		current_state_idx[1] += 1
		current_state_idx[0] += 1 
		current_state_idx[2] += 1 
		current_state = MATCH
		current_state_trans = match_list[0][1][:]
		path_string += 'M0->'
	else:
		# if start is insert
		path_likelihood = insert_list[0][0][nucDict[string[0]]]
		current_state_idx[1] += 1
		current_state_idx[0] += 1 
		current_state_idx[2] += 1 
		current_state = INSERT
		current_state_trans = insert_list[0][1][:]
		path_string += 'I0->'
	# now, all states should be possible
	for i in range(1, len(string)-1): # skips the \n at the end
		# decide which next state transition with the next emission is is most likely
		# TODO  ERRORS while looking ahead to far (easy fix)
		state_likelyhood_set = [ match_list[current_state_idx[0]+1][0][nucDict[string[i]]] * current_state_trans[0], insert_list[current_state_idx[1]+1][0][nucDict[string[i]]] * current_state_trans[1], current_state_trans[2] ]
		if current_state == INSERT:
			if ( insert_list[current_state_idx[1]+1][0][nucDict[string[i]]] * current_state_trans[1] == max(state_likelyhood_set) ):
				# we stay in the insert state we are in currently
				pass
			elif ( insert_list[current_state_idx[1]+1][0][nucDict[string[i]]] * current_state_trans[1] == max(state_likelyhood_set) ):
				# we transition into an match state
				current_state = MATCH
				current_state_idx[0] += 1
				current_state_trans = match_list[current_state_idx[0]][1][:]
				path_string += 'M' + str(i) + '->'
			else:
				# we transition into an delete state
				current_state = DELETE
				current_state_idx[2] += 1
				current_state_trans = delete_list[current_state_idx[2]][1][:]
				path_string += 'D' + str(i) + '->'
				i -= 1
		elif current_state == MATCH:
			if ( insert_list[current_state_idx[1]+1][0][nucDict[string[i]]] * current_state_trans[1] == max(state_likelyhood_set) ):
				# we transition into an insert state
				current_state_ = INSERT
				current_state_idx[1] += 1
				current_state_trans = insert_list[current_state_idx[1]][1][:]
				path_string += 'I' + str(i) + '->'
			elif ( insert_list[current_state_idx[1]+1][0][nucDict[string[i]]] * current_state_trans[1] == max(state_likelyhood_set) ):
				# we transition into an match state
				current_state = MATCH
				current_state_idx[0] += 1
				current_state_trans = match_list[current_state_idx[0]][1][:]
				path_string += 'M' + str(i) + '->'
			else:
				# we transition into an delete state
				current_state = DELETE
				current_state_idx[2] += 1
				current_state_trans = delete_list[current_state_idx[2]][1][:]
				path_string += 'D' + str(i) + '->'
				i -= 1
		elif current_state == DELETE:
			if ( insert_list[current_state_idx[1]+1][0][nucDict[string[i]]] * current_state_trans[1] == max(state_likelyhood_set) ):
				# we transition into an insert state
				current_state_ = INSERT
				current_state_idx[1] += 1
				current_state_trans = insert_list[current_state_idx[1]][1][:]
				path_string += 'I' + str(i) + '->'
			elif ( insert_list[current_state_idx[1]+1][0][nucDict[string[i]]] * current_state_trans[1] == max(state_likelyhood_set) ):
				# we transition into an match state
				current_state = MATCH
				current_state_idx[0] += 1
				current_state_trans = match_list[current_state_idx[0]][1][:]
				path_string += 'M' + str(i) + '->'
			else:
				# we transition into an delete state
				current_state = DELETE
				current_state_idx[2] += 1
				current_state_trans = delete_list[current_state_idx[2]][1][:]
				path_string += 'D' + str(i) + '->'
				i -= 1


					
	#print path_string




# ++++ main +++++++++++++++++++++++++++++++++++++++++++++++++++
def main():
	insert_list.append((pce, pct)) 												# add initial insert state
	for i in range(0, train_line_lenght): 										# iterate over chars in one column
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





	#for i, o in enumerate(match_list):
		#print 'M', i, ':', o[0], '=', sum(o[0]), o[1], '=', sum(o[1])


	# use viterbi on test files
	test_file = open(test_file_path)
	test_list = test_file.readlines()
	test_line_cnt = len(test_list)
	for i in range(1, test_line_cnt, 2):
		res = viterbi(test_list[i])
		print res






# =========================
if __name__ == "__main__":
	main()