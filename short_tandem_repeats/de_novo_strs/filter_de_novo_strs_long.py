#!/bin/python3

# filters de novo STRs with posterior > 0.05
# takes in *.all_mutations.tab created by MonSTR

import sys


infile = sys.argv[1]



fin = open(infile, 'r') 


wrote_header=False
for line in fin:

	if line.startswith('chrom'):
		if wrote_header:
			continue
		else:
			sys.stdout.write(line)
			wrote_header = True
			continue


	sline = line.split('\t')
	posterior = float(sline[7])

	if posterior < 0:
		sys.stdout.write(line)


fin.close()


