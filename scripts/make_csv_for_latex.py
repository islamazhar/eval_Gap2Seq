#count nr of gaps, total sun of N's in gaps.

import argparse
import sys
import csv
from itertools import izip
from itertools import izip_longest

#from itertools import *

# given an iterable of pairs return the key corresponding to the greatest value
def argmax(pairs):
    return max(pairs, key=lambda x:x[1])[0]
def argmin(pairs):
    return min(pairs, key=lambda x:x[1])[0]

def argmax_index(values):
    return argmax(enumerate(values))

def argmin_index(values):
    return argmin(enumerate(values))

def grouper(n, iterable, fillvalue=None):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

# def get_tool_difference(list_of_results):
# 	# Ceck so that it is numbers we are comparing
# 	try: 
# 		int(list_of_results[1])
# 	except:
# 		return list_of_results

# 	# only do difference for misassemblies
# 	if list_of_results[0] not in ['mismatches$^*$','short\_indels$^*$','long\_indels$^*$','misassemblies$^*$','local\_misassemblies$^*$' ]:
# 		return list_of_results

# 	# first two values are description and original assembly value
# 	new_list = [list_of_results[0],list_of_results[1]]
# 	print list_of_results
# 	# take difference
# 	original_value = int(list_of_results[1]) 
# 	for val in list_of_results[2:]:
# 		value = int(val)
# 		if value - original_value  < 0:
# 			new_list.append(str(value - original_value ))
# 		elif value -original_value > 0:
# 			new_list.append('+'+str(value - original_value ))
# 		else:
# 			new_list.append(value - original_value )
# 	return new_list

def make_bold(row,index):
	row[index] = '\\textbf{'+row[index]+'}'
	return row

def delete_tool_name(row):
	#print row
	return row[1:]

def hline(row):
	row[-1] = row[-1] + '\\hline '
	return row

def new_line(row):
	for i in range(len(row)):
		new_string = row[i] + '\\\ '
		print new_string
		row[i] = new_string
	print row
	return row

def sum_up_entrys_for_consensus(sum_stats, chunk_stats):
	for i in range(len(chunk_stats)):
		for j in range(len(chunk_stats[i])):
			if i == 0 or j <= 1:
				continue
			else:
				print 'i={0},j={1}'.format(i,j), int(chunk_stats[i][j])
				sum_stats[i][j-2] += int(chunk_stats[i][j]) 

def insert_sum(row,index):
	row.insert(index, '\\textbf{SUM}')
	return row

def bf_best_value(chunk_sum):
	print chunk_sum
	for i in range(len(chunk_sum[0])):
		stats = map(lambda x: x[i], chunk_sum[1:])
		if i == 0:
			continue
		elif i == 4:
			#print stats
			argmax = argmax_index(stats)
			#print argmax
			chunk_sum[argmax +1 ][i] = '\\textbf{'+str(chunk_sum[argmax + 1 ][i])+'}'
		else:
			argmin = argmin_index(stats)
			#print argmax
			chunk_sum[argmin +1 ][i] = '\\textbf{'+str(chunk_sum[argmin + 1 ][i])+'}'			

	return chunk_sum


def parse_quality(args):

	list_of_tuples_iter = list(csv.reader(open(args.quality,'rb')))

	#header is special case
	header = list_of_tuples_iter[0] #map(lambda x: [x],list_of_tuples_iter[0])
	#header = [string.replace('_','\_')for string in header]
	#header = [string.replace('#','')for string in header]
	a = izip(*header)

	# the rest of the rows with the tool results
	csv_file = open(args.outprefix, "wb")

	# initialize the total sum of stats vector
	chunk_sum = []
	for k in range(6):
		chunk_sum.append([0,0,0,0,0,0])

	for i, tool_chunk in enumerate(grouper(5, list_of_tuples_iter[1:])):
	#for tool_chunk in csv.reader(open(args.quality,'rb')):
		tool_chunk = list(tool_chunk)
		tool_chunk.insert(0,header)
		print tool_chunk

		# remowe rhodo entry for S-aurens
		if tool_chunk[-1][-1] =='-':
			continue

		if i == 0:
			tool_chunk[4][0] = tool_chunk[4][0].replace('_','-') 
			tool_chunk[5][0] = tool_chunk[5][0].replace('_','-')	
			tool_chunk = [make_bold(x,1) for x in tool_chunk]
			chunk_sum[0] = header[2:]
			sum_up_entrys_for_consensus(chunk_sum, tool_chunk)

		else:
			sum_up_entrys_for_consensus(chunk_sum, tool_chunk)
			tool_chunk = map(delete_tool_name, tool_chunk)
			tool_chunk = [make_bold(x,0) for x in tool_chunk]


		#print tool_chunk[-1]
		tool_chunk[-1] = new_line(tool_chunk[-1])
		#print tool_chunk[-1]
		tool_chunk[-1] = hline(tool_chunk[-1])

		a = izip(* tool_chunk)

		# for i in a:
		# 	print i
		csv.writer(csv_file, delimiter='&').writerows(a)

	print 'sum'
	print chunk_sum
	chunk_sum = [insert_sum(x,0) for x in chunk_sum]
	bf_best_value(chunk_sum)
	chunk_sum[-1] = new_line(map(lambda x: str(x), chunk_sum[-1]))
	a = izip(* chunk_sum)
	csv.writer(csv_file, delimiter='&').writerows(a)




def parse_performance(args):

	list_of_tuples_iter = list(csv.reader(open(args.performance,'rb')))

	# the rest of the rows with the tool results
	csv_file = open(args.outprefix, "wb")
	for i, tool_chunk in enumerate(grouper(5, list_of_tuples_iter[1:])):
		
		tool_chunk = list(tool_chunk)
		tool_chunk.insert(0,['tool','assembly','User time','Wall clock time','memory'])

		# remowe rhodo entry for S-aurens
		if tool_chunk[-1][-1] =='-':
			continue

		if i == 0:
			tool_chunk[4][0] = tool_chunk[4][0].replace('_','-') 
			tool_chunk[5][0] = tool_chunk[5][0].replace('_','-')	
			tool_chunk = [make_bold(x,1) for x in tool_chunk]

		else:
			tool_chunk = map(delete_tool_name, tool_chunk)
			tool_chunk = [make_bold(x,0) for x in tool_chunk]

		#print tool_chunk[-1]
		tool_chunk[-1] = new_line(tool_chunk[-1])
		#print tool_chunk[-1]
		tool_chunk[-1] = hline(tool_chunk[-1])

		a = izip(* tool_chunk)
		csv.writer(csv_file, delimiter='&').writerows(a)

def main(args):
	if args.quality:
		qual_dict = parse_quality(args)
	elif args.performance:
		perf_dict = parse_performance(args)
	


if __name__ == '__main__':
    ##
    # Take care of input
	parser = argparse.ArgumentParser(description="Calculate number of gaps and total number of N's in a scaffold.fasta file")
	parser.add_argument('-q', dest='quality',type=str, help='A quality csv file.')
	parser.add_argument('-p', dest='performance', type=str, help='A performance csv file.')

	parser.add_argument('-o', dest='outprefix', type=str, help='prefix to outfile. (Outfile is selected fields formatted for latex table')

	args = parser.parse_args()
	main(args)