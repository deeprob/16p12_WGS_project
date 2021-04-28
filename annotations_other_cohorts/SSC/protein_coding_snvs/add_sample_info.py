#!/bin/python3


import gzip
import sys
import vcf

vcf_with_annotations_filename = sys.argv[1]
vcf_with_samples_filename = sys.argv[2]
error_filename = sys.argv[3] # will write error messages here

def get_record_at(chrom, pos, alt_allele, ref_allele):
	for rec in vcfr.fetch(chrom, int(pos)-1, int(pos)+1):
		alt_alleles = rec.ALT
		
		# if start positions don't match go to next one
		if int(pos) != rec.POS:
			continue
		
		# if ref alleles don't match go to next one
		if ref_allele != rec.REF:
			continue
		
		# check that alt allele in anno record
		for a in alt_alleles:
			if a == alt_allele:
				return rec
		
	# will need to add these in manually later
	with open(error_filename, 'a') as f:
		error_message = 'no record at {} {} {} {}\n'.format(chrom, pos, ref_allele, alt_allele)
		f.write(error_message) 

def get_alt_allele_int(rec, alt_allele):
	sample_alt_alleles = rec.ALT
	for i in range(len(sample_alt_alleles)):
		a = sample_alt_alleles[i]
		if a == alt_allele:
			return i + 1

def quality_score(rec, data, alt_allele_int):
	# returns True if PASS and False if FAIL
	# check if DP is none
	if data.DP is None:
		return False
	# check if FORMAT/DP>=8
	if data.DP < 8:
		return False
	# check if FORMAT/AD[:1]>0
	if data.AD[alt_allele_int] == 0:
		return False
	# check if (FORMAT/AD[:1])/(FMT/DP)>=0.25
	if data.AD[alt_allele_int] / data.DP < 0.25:
		return False
	# check if (FMT/AD[:1])/(FMT/DP)<=0.75
	if data.AD[alt_allele_int] / data.DP > 0.75:
		return False
	# check if QUAL/(FMT/AD[:1])>=1.5
	if rec.QUAL / data.AD[alt_allele_int] < 1.5:
		return False
	return True

# open error message filename just to empty it out
with open(error_filename, 'w') as f:
	pass

# load vcf with samples
vcfr = vcf.Reader(filename = vcf_with_samples_filename)

# load vcf with annotations
vcfanno = gzip.open(vcf_with_annotations_filename, 'r')

for line in vcfanno:
	line = line.decode()
		
	# header
	if line.startswith('#'):
		continue
		
	# split by tab
	sline = line.strip().split('\t')
	
	chrom      = sline[0]
	pos        = sline[1]
	ref_allele = sline[3]
	alt_allele = sline[4]
	
	rec = get_record_at(chrom, pos, alt_allele, ref_allele)
	 
	# if no record is found continue (put it in manually later)
	if rec is None:
		continue
	
	# sample VCF can have multi-allele records
	# find out which allele is the one we want
	alt_allele_int = get_alt_allele_int(rec, alt_allele)
	
	# for each sample with the alternate allele
	# check that the sample quality score passes
	# then write out
	sample_gts = rec.samples
	for sample_gt in sample_gts:
		sample = sample_gt.sample
		data   = sample_gt.data
		
		# check if has allele we want
		has_alt_allele = str(alt_allele_int) in data.GT
		
		# if sample doesn't have alt allele the skip
		if not has_alt_allele:
			continue
		
		# check if samples passes quality metrics
		quality = quality_score(rec, data, alt_allele_int)	
		
		# if sample doesn't pass quality metrics then skip
		if not quality:
			continue

		# construct outline and write out
		outline = '\t'.join(sline) + '\t' + sample + '\n'
		sys.stdout.write(outline)




	
vcfanno.close()




