import argparse
import os,sys
import re

def parse_quast_stdout(infile):
	pass

# def parse_quast_misassemblies(infile, out_string):
# 	lines = map( lambda line: line.strip().split('  '), infile)
# 	filtered_lines=[]
# 	for i,line in enumerate(lines):
# 		content = filter(lambda x: x != '',line )
# 		rem_whitespace_content=[]
# 		if i in [10,12,13]:
# 			for col in content:
# 				rem_whitespace_content.append(col.replace(" ", ""))
# 			out_string +=rem_whitespace_content[2]+','
# 	return out_string

def parse_quast_report(infile, out_string):
	lines = map( lambda line: line.strip().split('  '), infile)
	filtered_lines=[]
	for i,line in enumerate(lines):
		content = filter(lambda x: x != '',line )
		rem_whitespace_content=[]
		if i in [21,26,34]:
			for col in content:
				rem_whitespace_content.append(col.replace(" ", ""))
			if i == 21:
				out_string += rem_whitespace_content[2]
			else:
				out_string += ','+rem_whitespace_content[2]
	return out_string


def parse_gnu_time(infile):
	performance_string = ''
	for l in infile:
		#print l
		usertime =  re.search('User time \(seconds\): [\d.]+', l)
		wct= re.search('Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): [\d.:]+', l) 
		mem = re.search('Maximum resident set size \(kbytes\): [\d.:]+', l) 
		if usertime:
			performance_string += usertime.group().split(':')[1].strip()
		if wct:
			performance_string += ','+wct.group().split()[7]
		if mem:
			mem_tmp = int(mem.group().split()[5])
			memory = mem_tmp / 4 
			performance_string += ','+str(memory)

	return performance_string

def main(args):
	if not os.path.exists(args.outpath):
		os.mkdir(args.outpath)
	#parse_quast_stdout(open(args.quast_stdout,'r'))
	quality_string = ''
	try:
		#quality_string = parse_quast_misassemblies(open(args.quast_misassemblies,'r'),quality_string)
		quality_string = parse_quast_report(open(args.quast_report,'r'), quality_string)
		print >> open(os.path.join(args.outpath,'gap_filling_quality_eval.csv'),'w' ), quality_string
	except IOError:
		print '{0} does not exist. skipping...'.format(args.quast_misassemblies)
		print >> open(os.path.join(args.outpath,'gap_filling_quality_eval.csv'),'w' ), '-,-,-,-,-,-'

	# quality_string = parse_quast_misassemblies(open(args.quast_misassemblies,'r'),quality_string)
	# quality_string = parse_quast_report(open(args.quast_report,'r'), quality_string)
	# print >> open(os.path.join(args.outpath,'gap_filling_quality_eval.csv'),'w' ), quality_string
	try:
		performance_string = parse_gnu_time(open(args.gnu_time,'r'))
		print >> open(os.path.join(args.outpath,'gap_filling_performance_eval.csv'),'w' ), performance_string
	except IOError:
		print '{0} does not exist. skipping...'.format(args.gnu_time)
		print >> open(os.path.join(args.outpath,'gap_filling_performance_eval.csv'),'w' ), '-,-,-'
	


if __name__ == '__main__':
    ##
    # Take care of input
	parser = argparse.ArgumentParser(description="Evaluate gaps filled by the Salmela-Ford algorithm.")
	parser.add_argument('quast_stdout', type=str, help='A quast A quast contig_reports/contig_reports.stdout file.')
	parser.add_argument('quast_misassemblies', type=str, help='A quast contig_reports/missassemblies.txt file')
	parser.add_argument('quast_report', type=str, help='A quast report.txt file.')
	parser.add_argument('gnu_time', type=str, help='An outsahlin.stderr file with gnu time.')

	parser.add_argument('outpath', type=str, help='Folder for output.')
	# parser.add_argument('sigma', type=float, help='Stddev of library')
	# parser.add_argument('cov', type=float, help='Mean coverage of library')
	# parser.add_argument('r', type=float, help='read length  of library')
	# parser.add_argument('s', type=float, help='Maximum allowed softclipped bases')

	args = parser.parse_args()
	main(args)