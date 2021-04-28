#!/bin/python3


import sys
import gzip


infile = 'vcfs/genes/chr22.hg38_multianno.vcf.gz'


infile = sys.argv[1]


def get_info_field(info_line, info_field):
        # returns the info field
        sinfo = info_line.split(';')
        for item in sinfo:
                if item.startswith(info_field + '='):
                        length_info_field = len(info_field)
                        return item[length_info_field + 1:]
        return '.'

def get_exonic_function(info_line):
        # returns the exonic function
        exonic_function = get_info_field(info, 'ExonicFunc.refGene')
        # categorize into large groups
        if exonic_function in missense_variant_functions:
        	return 'missense'
        if exonic_function in lof_variant_functions:
        	return 'lof'
        if exonic_function in frameshift_variant_function:
        	return 'frameshift'
        return exonic_function

f = gzip.open(infile)


new_header_line = '##INFO=<ID=variant_function,Number=.,Type=String,Description="Protein coding variant function">\n'


locations_keep = ['exonic', 'splicing', 'exonic\\x3bsplicing']
missense_variant_functions = ['nonsynonymous_SNV', 'nonframeshift_insertion',
       'nonframeshift_deletion']
lof_variant_functions = ['stopgain', 'stoploss']
frameshift_variant_function = ['frameshift_insertion', 'frameshift_deletion']
variant_functions_to_keep = ['missense', 'lof', 'frameshift', 'splice']

for line in f:
	line = line.decode()

	# write new header line
	if line.startswith('#CHROM'):
		sys.stdout.write(new_header_line)
		sys.stdout.write(line)
		continue

	if line.startswith('#'):
		sys.stdout.write(line)
		continue

	# split by tab
	sline = line.strip().split('\t')
	qual  = sline[5]
	info  = sline[7]

	# if qual is < 50 then skip
	if float(qual) < 50:
		continue

	# get gnomad frequency
	freq = get_info_field(info, 'gnomad_AF')
	freq = float(freq)

	# if frequency > 0.001 then skip
	if freq > 0.001:
		continue

	# get location
	location = get_info_field(info, 'Func.refGene')
	if location not in locations_keep:
		continue

	# get variant function
	if location == 'splicing':
		exonic_function = 'splice'
	else:
		exonic_function = get_exonic_function(info)

	# if exonic function is synonymous then skip
	if exonic_function not in variant_functions_to_keep:
		continue

	# add variant function to line
	outline = line[:-1] + ';variant_function=' + exonic_function + '\n'

	sys.stdout.write(outline)


f.close()







