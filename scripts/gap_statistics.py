#count nr of gaps, total sun of N's in gaps.

import re
import argparse
import sys

def fasta_iter(fasta_file):
    """
        Reads a fasta file into memory.

        Arguments:
        fasta_file - A python file object. The file should be in 
        fasta format.

        Returns:
            an iterator over accession, sequence.

    """  

    k = 0
    temp = []
    accession = ''
    for line in fasta_file:
        if line[0] == '>' and k == 0:
            accession = line[1:].strip().split()[0]
            k += 1
        elif line[0] == '>':
            temp = ''.join(temp)
            yield accession, temp
            temp = []
            accession = line[1:].strip().split()[0]
        else:
            temp.append(line.strip())
    
    temp = ''.join(temp)
    yield accession, temp

def main(args):
	try:
		scf_fa_file = open(args.infile,'r')
	except IOError:
		sys.stderr.write('{0} does not exist. skipping...\n'.format(args.infile))
		print '{0},{1}'.format('-','-')
		sys.exit()

	nr_gaps = 0
	tot_gap_length = 0
	gap_distr = []
	for acc, seq in fasta_iter(scf_fa_file):
		gap_list = re.findall('[Nn]+', seq)
		gap_lengths = map(lambda x: len(x), gap_list)
		gap_distr += gap_lengths
		tot_gap_length += sum(gap_lengths)
		nr_gaps += len(gap_lengths)

	#print 'nr_gaps\ttot_gap_length'
	print '{0},{1}'.format(nr_gaps,tot_gap_length)
	#print gap_distr

if __name__ == '__main__':
    ##
    # Take care of input
	parser = argparse.ArgumentParser(description="Calculate number of gaps and total number of N's in a scaffold.fasta file")
	parser.add_argument('infile', type=str, help='A scaffold file in fasta format')
	# parser.add_argument('sigma', type=float, help='Stddev of library')
	# parser.add_argument('cov', type=float, help='Mean coverage of library')
	# parser.add_argument('r', type=float, help='read length  of library')
	# parser.add_argument('s', type=float, help='Maximum allowed softclipped bases')

	args = parser.parse_args()
	main(args)