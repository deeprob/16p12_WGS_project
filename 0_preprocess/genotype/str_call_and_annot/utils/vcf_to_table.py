#!/bin/python3

import gzip
import pandas as pd
import sys

bedfilename = '/data5/16p12_WGS/structural_variants/vcf_callers/hg19_ver13_1.bed'
vcf_filename='output/16p12_cohort.vcf.gz'


bed = pd.read_csv(bedfilename, sep='\t', header=None)
bed['variant_id'] = bed[0] + '_' + bed[1].astype(str)
bed = bed.set_index('variant_id')

fin = gzip.open(vcf_filename)
for line in fin:
	line = line.decode()
	
	# write header
	if line.startswith('#CHROM'):
		sline = line.strip().split('\t')
		samples = sline[9:]
		outline = 'chrom\tpos\t'
		outline = outline + '\t'.join(samples) + '\n'
		sys.stdout.write(outline)
		continue
	
	if line.startswith('#'):
		continue
	
	sline = line.strip().split('\t')

	chrom = sline[0]
	pos = sline[1]
	ref_allele = sline[3]
	alt_alleles = sline[4]
	gt_infos = sline[9:]
	
	# get period of str
	variant_id = chrom + '_' + str(pos)
	period = bed.loc[variant_id, 3]
	
	# get lengths of all of the alleles
	if alt_alleles == '.':
		alleles = [ref_allele]
	else:
		alleles = [ref_allele] + alt_alleles.split(',')
	alleles_lengths = [len(s)/period for s in alleles]
	
	# strip extra info from gts
	gt_infos = [s.split(':')[0] for s in gt_infos]
	
	# process gts and create outline
	outline = chrom + '\t' + str(pos) + '\t'
	for gt_info in gt_infos:
		if gt_info == '.':
			outline = outline + '.\t'
			continue
		gt = gt_info.split('/')
		lengths = [alleles_lengths[int(s)] for s in gt]
		lengths = [str(int(s)) for s in lengths]
		outline = outline + ','.join(lengths) + '\t'
	
	# remove trailing tab and add new line
	outline = outline[:-1] + '\n'
	sys.stdout.write(outline)



fin.close()

