#!/bin/python3


import sys
import pandas as pd

filename = sys.argv[1]
pedfilename = sys.argv[2]

# load in sex to check X chromosome
ped = pd.read_csv(pedfilename, sep='\t', header=None)
ped = ped.set_index(1, drop=False)

fin = open(filename, 'r')

for line in fin:

	# header 
	if line.startswith('chrom'):
		sys.stdout.write(line)
		continue


	sline = line.split('\t')
	
	chrom = sline[0]
	child_gts = sline[16]
	child_gts = child_gts.split(',')
	child = sline[5]

	child_x_ploidy = ped.loc[child, 4]

	# if child is male and mutation is on X chromosome
	# keep even though it's homoxygous
	if child_x_ploidy == 1 and chrom == 'X':
		sys.stdout.write(line)
		continue

	if child_gts[0] != child_gts[1]:
		sys.stdout.write(line)

fin.close()




