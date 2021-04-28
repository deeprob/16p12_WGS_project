#!/bin/python3


# filters out variants with QUAL < 50 or gnomad > 0.001
# and only keep variants that are exonic or splicing

import gzip
import sys

invcf = 'SPARK_pilot1379ind.ColumbiaJointCall.annovar.hg19_multianno.vcf.gz'

invcf = sys.argv[1]



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
	if dp == '.':
		return False
	dp = int(dp)
	ad = int(ad.split(',')[1])
	qual = float(qual)
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



locations_keep = ['exonic', 'splicing', 'exonic-splicing']
missense_variant_functions = ['nonsynonymous_SNV', 'nonframeshift_insertion',
       'nonframeshift_deletion']
lof_variant_functions = ['stopgain', 'stoploss']
frameshift_variant_function = ['frameshift_insertion', 'frameshift_deletion']
variant_functions_to_keep = ['missense', 'lof', 'frameshift', 'splice']


f = gzip.open(invcf, 'r')
for line in f:
	line = line.decode()
	#
	# get sample names
	if line.startswith('#CHROM'):
		sline = line.strip().split('\t')
		samples = sline[9:]
		# sys.stdout.write(new_header_line)
		# sys.stdout.write(line)
		continue
	#
	# skip header
	if line.startswith('#'):
		# sys.stdout.write(line)
		continue
	#
	# split by tab
	sline = line.strip().split('\t')
	#
	chrom      = sline[0]
	pos        = sline[1]
	ref_allele = sline[3]
	alt_allele = sline[4]
	qual       = sline[5]
	info       = sline[7]
	format     = sline[8]
	gt_infos   = sline[9:]
	#
	# if qual < 50 then skip
	if float(qual) < 50:
		continue
	#
	# get gnomad frequency
	gnomad_freq = get_info_field(info, 'gnomAD_genome_ALL')
	# 
	# convert gnomad freq from string to float
	if gnomad_freq == '.':
		gnomad_freq = 0
	else:
		gnomad_freq = float(gnomad_freq)
	#
	# if gnomad freq > 0.001 then skip
	if gnomad_freq > 0.001:
		continue
	#
	# get location
	location = get_info_field(info, 'Func.refGene')
	if location not in locations_keep:
		continue
	#
	# get variant function
	if location == 'splicing':
		exonic_function = 'splice'
	else:
		exonic_function = get_exonic_function(info)
	#
	# if exonic function is synonymous then skip
	if exonic_function not in variant_functions_to_keep:
		continue
	#
	# loop through each sample
	for i in range(len(samples)):
		sample  = samples[i]
		gt_info = gt_infos[i]
		#
		# get GT 
		gt = get_format_field(format, 'GT', gt_info)
		# 
		# if sample doesn't have alt allele then skip
		if '1' not in gt:
			continue
		# get DP and AD
		dp = get_format_field(format, 'DP', gt_info)
		ad = get_format_field(format, 'AD', gt_info)
		#
		# get quality
		quality = quality_score(qual, dp, ad)
		# 
		# if sample doesn't pass quality then skip
		if not quality:
			continue
		#
		# construct outline
		new_sline = sline[:9]
		#
		# add new info field
		new_info = info + ';variant_function=' + exonic_function
		new_sline[7] = new_info
		#
		# append sample info
		new_sline.append(gt_info)
		new_sline.append(sample)
		#
		# write out new line
		new_line = '\t'.join(new_sline) + '\n'
		sys.stdout.write(new_line)









