#count nr of gaps, total sun of N's in gaps.

import re
import argparse
import sys
import os


def main(args):
	experiment_out = os.path.join(args.infolder, 'experiment-{0}/gap_filling_quality_eval.csv'.format(1))
	try:
		result= open(experiment_out,'r')
	except IOError:
		sys.stderr.write('{0} does not exist. skipping...\n'.format(experiment_out))
		sys.exit()

	sum_result = [0,0,0,0]

	for i in range(1,args.nr_experiments+1):
		experiment_out = os.path.join(args.infolder, 'experiment-{0}/gap_filling_quality_eval.csv'.format(i))
		result = open(experiment_out,'r').readline().strip().split(',')
		res = map(lambda x: int(x), result)
		print res
		for j in range(len(res)):
			sum_result[j] += res[j]

	average = map(lambda x: int(round(x/float(args.nr_experiments),0)) , sum_result)
	outstring = ','.join([ str(x) for x in average])
	print outstring
	outfile = os.path.join(args.infolder, 'gap_filling_quality_eval.csv')
	print >> open(outfile, 'w' ), outstring


if __name__ == '__main__':
    ##
    # Take care of input
	parser = argparse.ArgumentParser(description="Calculate number of gaps and total number of N's in a scaffold.fasta file")
	parser.add_argument('infolder', type=str, help='folder for gap tool experiments')
	parser.add_argument('nr_experiments',type=int, help='How many experiments that was run.')


	args = parser.parse_args()
	main(args)