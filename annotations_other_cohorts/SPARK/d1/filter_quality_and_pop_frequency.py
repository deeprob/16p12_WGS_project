#!/bin/python3


import gzip
import sys
import vcf

invcf = sys.argv[1]


def get_format_field(format, field, gt_info):
	sformat = format.split(':')

	# get index of field
	for i in range(len(sformat)):
		item = sformat[i]
		if item == field:
			field_index = i

	return gt_info.split(':')[field_index]

def quality_score(qual, dp, ad):
	# returns True if PASS and False if FAIL

	dp = int(dp)
	ad = int(ad.split(',')[1])
	qual = int(qual)

	# check if FORMAT/DP>=8
	if dp < 8:
		return False
	# check if FORMAT/AD[:1]>0
	if ad == 0:
		return False
	# check if (FORMAT/AD[:1])/(FMT/DP)>=0.25
	if ad / dp < 0.25:
		return False
	# check if (FMT/AD[:1])/(FMT/DP)<=0.75
	if ad / dp > 0.75:
		return False
	# check if QUAL/(FMT/AD[:1])>=1.5
	if qual / ad < 1.5:
		return False
	return True

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

locations_keep = ['exonic', 'splicing', 'exonic\\x3bsplicing']
missense_variant_functions = ['nonsynonymous_SNV', 'nonframeshift_insertion',
       'nonframeshift_deletion']
lof_variant_functions = ['stopgain', 'stoploss']
frameshift_variant_function = ['frameshift_insertion', 'frameshift_deletion']
variant_functions_to_keep = ['missense', 'lof', 'frameshift', 'splice']

# load vcf
vcfanno = gzip.open(invcf, 'r')

for line in vcfanno:
	line = line.decode()

	if line.startswith('#CHROM'):
		sline = line.strip().split('\t')
		sample = sline[-1]

	# header
	if line.startswith('#'):
		continue
		
	# split by tab
	sline = line.strip().split('\t')
	
	chrom      = sline[0]
	pos        = sline[1]
	ref_allele = sline[3]
	alt_allele = sline[4]
	qual       = sline[5]
	info       = sline[7]
	format     = sline[8]
	gt_info    = sline[9]

	# if qual < 50 skip
	if float(qual) < 50:
		continue

	# check frequency of alternative allele
	freq = get_info_field(info, 'AF')

	# if frequency is greater than 0.001 skip
	if float(freq) > 0.001:
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

	# get DP and AD
	dp = get_format_field(format, 'DP', gt_info)
	ad = get_format_field(format, 'AD', gt_info)

	# check quality of record
	quality = quality_score(qual, dp, ad)

	# if record doesn't pass quality standrards then skip
	if not quality:
		continue

	# construct outline and write out
	new_sline    = sline
	new_info     = info + ';variant_function=' + exonic_function
	new_sline[7] = new_info
	outline      = '\t'.join(new_sline) + '\t' + sample + '\n'
	sys.stdout.write(outline)



	
vcfanno.close()






