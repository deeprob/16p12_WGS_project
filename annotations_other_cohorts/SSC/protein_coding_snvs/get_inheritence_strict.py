#!/bin/python3

# this script adds inheritence patterns

import sys
import pandas as pd
import vcf

infilename = 'vcfs/protein_coding/all_chromosomes.select_annotations.leouf.tsv'
sample_map_filename = '../documentation/nygc_sfari_id_map.csv'
outfilename = 'vcfs/protein_coding/all_chromosomes.select_annotations.leouf.new_inh.tsv'

# load in samp mapp
sampmapp = pd.read_csv(sample_map_filename)
sampmapp = sampmapp.set_index('SFARI ID')['Repository Id'].to_dict()

# load in variant
df = pd.read_csv(infilename, sep='\t', low_memory=False)

df['variant_id'] = df.Chr + '_' + df.Pos.astype(str) + '_' + df.Alt
df['inheritence_strict'] = '.'
df['Family'] = df['Family'].astype(str)


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



current_chrom = ''
current_var_id = ''
for i, row in df.iterrows():
	sample     = row['Sample']
	family     = row['Family']
	chrom      = row['Chr']
	pos        = row['Pos']
	ref_allele = row['Ref']
	alt_allele = row['Alt']
	var_id     = row['variant_id']

	# if sample is a parent then skip
	if sample.endswith('mo') or sample.endswith('fa'):
		continue

	# get mother and father
	mother = sampmapp[family + '.mo']
	father = sampmapp[family + '.fa']

	# if new chrom then load in correct vcf file
	if chrom != current_chrom:
		print(chrom)
		current_chrom = chrom
		vcf_filename = 'vcfs/norm/{}.vcf.gz'.format(chrom)
		vcfr = vcf.Reader(filename = vcf_filename)

	# if mother or father are not in vcf then skip
	if (mother not in vcfr.samples) or (father not in vcfr.samples):
		continue

	# if new location then get vcf record
	if var_id != current_var_id:
		# print(var_id)
		current_var_id = var_id
		record = get_record_at(chrom, pos, alt_allele, ref_allele)

	# safety check
	if len(record.ALT) != 1:
		sys.stderr.write(record, ' more than 1 alt allele')
		break

	# get call for mother and father
	mother_call = record.genotype(mother).data
	father_call = record.genotype(father).data

	# get quality for mother and father
	mother_quality = quality_score(record, mother_call)
	father_quality = quality_score(record, father_call)

	# if either mother or the father do not pass quality then skip
	if (not mother_quality) or (not father_quality):
		continue

	# finally get inheritence
	if '1' in mother_call.GT and '1' in father_call.GT:
		inheritence = 'both'
	elif '1' in mother_call.GT:
		inheritence = 'mother'
	elif '1' in father_call.GT:
		inheritence = 'father'
	else:
		inheritence = 'de novo'

	df.at[i, 'inheritence_strict'] = inheritence



# save to file
df.to_csv(outfilename, sep='\t', index=False)









