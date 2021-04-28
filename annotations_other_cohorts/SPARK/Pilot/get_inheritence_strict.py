#!/bin/python3

# this script adds inheritence patterns

import sys
import vcf
import pandas as pd

infilename = 'SPARK_pilot1379ind.ColumbiaJointCall.filtered.select_annotations.txt'
sample_map_filename = '/data4/SPARK/Pilot/documentation/1379ind_masterTable_20181016.txt'
vcf_filename = 'SPARK_pilot1379ind.ColumbiaJointCall.norm.vcf.gz'



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
		if alt_allele in alt_alleles:
			return rec

def quality_score(rec, data):
	# returns True if PASS and False if FAIL
	# check if GT is empty
	if data.GT == './.':
		return False
	# check if DP is none
	if data.DP is None:
		return False
	# check if FORMAT/DP>=8
	if data.DP < 8:
		return False
	# if Reference allele no need to check further
	if data.GT == '0/0':
		return True
	# check if FORMAT/AD[:1]>0
	if data.AD[1] == 0:
		return False
	# check if (FORMAT/AD[:1])/(FMT/DP)>=0.25
	if data.AD[1] / data.DP < 0.25:
		return False
	# check if (FMT/AD[:1])/(FMT/DP)<=0.75
	if data.AD[1] / data.DP > 0.75:
		return False
	# check if QUAL/(FMT/AD[:1])>=1.5
	if rec.QUAL / data.AD[1] < 1.5:
		return False
	return True

# open vcf file
vcfr = vcf.Reader(filename = vcf_filename)

# load in samp mapp
names = ['SPID','SFID','Father','Mother','Gender','Phenotype','IDAT','Array_FinalReports','Bam','Comment']
sampmapp = pd.read_csv(sample_map_filename, sep='\t', names=names, skiprows=1)
sampmapp = sampmapp.set_index('SPID', drop=False)

# get list of all fathers and all mothers
fathers = list(sampmapp['Father'].unique())
mothers = list(sampmapp['Mother'].unique())

f = open(infilename, 'r')

current_var_id = ''

for line in f:
	#
	# header
	if line.startswith('Family'):
		new_line = line.strip()
		new_line = new_line + '\tinheritence_strict\n'
		sys.stdout.write(new_line)
		continue
	#
	sline = line.strip().split('\t')
	#
	family      = sline[0]
	sample      = sline[1]
	chrom       = sline[2]
	pos         = sline[3]
	ref_allele  = sline[4]
	alt_allele  = sline[5]
	#
	# if sample is a parent then skip
	if sample in fathers or sample in mothers:
		new_line = line.strip()
		new_line = new_line + '\t.\n'
		sys.stdout.write(new_line)
		continue
	#
	# find out who is father and who is mother
	mother = sampmapp.loc[sample, 'Mother']
	father = sampmapp.loc[sample, 'Father']
	#
	# if mother or father are not in vcf then skip
	if (mother not in vcfr.samples) or (father not in vcfr.samples):
		new_line = line.strip()
		new_line = new_line + '\t.\n'
		sys.stdout.write(new_line)
		continue
	#
	# if new location then get vcf record
	var_id = chrom + '_' + pos + '_' + alt_allele
	if var_id != current_var_id:
		# print(var_id)
		current_var_id = var_id
		record = get_record_at(chrom, pos, alt_allele, ref_allele)
	#
	# safety check
	if len(record.ALT) != 1:
		sys.stderr.write(record, ' more than 1 alt allele')
		break
	#
	# get call for mother and father
	mother_call = record.genotype(mother).data
	father_call = record.genotype(father).data
	#
	# get quality for mother and father
	mother_quality = quality_score(record, mother_call)
	father_quality = quality_score(record, father_call)
	#
	# if either mother or the father do not pass quality then skip
	if (not mother_quality) or (not father_quality):
		new_line = line.strip()
		new_line = new_line + '\t.\n'
		sys.stdout.write(new_line)
		continue
	#
	# finally get inheritence
	if '1' in mother_call.GT and '1' in father_call.GT:
		inheritence = 'both'
	elif '1' in mother_call.GT:
		inheritence = 'mother'
	elif '1' in father_call.GT:
		inheritence = 'father'
	else:
		inheritence = 'de novo'
	#
	# construct new line
	new_line = line.strip()
	new_line = new_line + '\t' + inheritence + '\n'
	sys.stdout.write(new_line)


f.close()










