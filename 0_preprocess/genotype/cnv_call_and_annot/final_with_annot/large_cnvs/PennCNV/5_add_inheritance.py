import pandas as pd
import subprocess
import datetime

#####################################################################################################
# IMPORTANT: This script does not work on qingyu!!
# bcftools error:
# bcftools: /lib/x86_64-linux-gnu/libm.so.6: version `GLIBC_2.29' not found (required by bcftools)
#####################################################################################################

# Step to avoid truncation of long columns - relevant for getting filenames
pd.options.display.max_colwidth = 999

# Use the lookup table to find samples that have the same CNV (defined as 50% reciprocal overlap)
calls = pd.read_csv('call_tables/3_annotate_gencode.txt', sep = '\t')
lookup = pd.read_csv('bed_files/4_inheritance_lookup.bed', sep = '\t',
	names = ['Chr', 'Start', 'End', 'Type', 'Sample', 'variant_id', 'Chr2', 'Start2', 'End2', 'Type2', 'Sample2', 'variant_id2'])

fam = pd.read_csv('../../2022_02_14/Feb_2022.fam', sep = '\t', names = ['Family', 'Sample', 'Father', 'Mother', 'Sex', 'Carrier_Status'])
penncnv = pd.read_csv('../../2022_02_14/PennCNV_calls_combined_final.csv')

# Find 50% reciprocal overlap calls in the final output and return the samples they belong to
def overlap_samps(vid):
	overlaps = lookup[(lookup.variant_id==vid) & (lookup.variant_id2!=vid)]
	samps = list(overlaps.Sample2.unique())

	return samps

# Function to check for reciprocal overlap
def recip_overlap(start1, end1, start2, end2, length1, length2, threshold = 0.5):
	overlap = (min(end1, end2) - max(start1, start2))
	if overlap >= (threshold*max(length1, length2)):
		return True
	return False

# Look through the raw PennCNV output
def penncnv_lookup(sample, chr, start, end, length, svtype):
	# Get lines for parent from PennCNV calls
	df = penncnv[penncnv.PatientID==sample]
	if df.shape[0]==0:
		return False
	# Get correct chromosome and type
	# Change type because PennCNV calls have type differently
	new_type = svtype.replace('<', "").replace('>', "").lower()
	df = df[(df.Chromosome==chr) & (df.Type==new_type)]
	if df.shape[0]==0:
		return False
	# Get correct location
	search_start = start-length
	search_end = start+(length/2)
	df = df[(df.Start >= search_start) & (df.Start <= search_end)]
	if df.shape[0]==0:
		return False
	# Check for 50% reciprocal overlap
	df = df[df.apply(lambda row: recip_overlap(start, end, row['Start'], row['End'], length, row['Length']), axis = 1)]
	if df.shape[0]>0:
		return True
	return False

def inheritance(row, father = False, mother = False):
	iid = row['PatientID']

	# Get family information from FAM file
	mother_code = fam[fam.Sample==iid]['Mother'].to_string(index = False, header = False).strip()
	father_code = fam[fam.Sample==iid]['Father'].to_string(index = False, header = False).strip()

	# Check if we have CNV calls for the parents
	if mother_code == '0' and father_code == '0':
		return '.'
	if mother_code not in penncnv.PatientID.to_list() and father_code not in penncnv.PatientID.to_list():
		return '.'

	# Get carrier status of parents
	if mother_code == '0':
		mother_return = '0'
	elif fam[fam.Sample==mother_code]['Carrier_Status'].to_string(index = False).strip()=='unknown':
		mother_return = 'M'
	elif int(fam[fam.Sample==mother_code]['Carrier_Status'])==1:
		mother_return = 'MC'
	elif int(fam[fam.Sample==mother_code]['Carrier_Status'])==0:
		mother_return = 'MNC'
	else:
		mother_return = 'M'

	if father_code == '0':
		father_return = '0'
	elif fam[fam.Sample==father_code]['Carrier_Status'].to_string(index = False).strip()=='unknown':
		father_return = 'F'
	elif int(fam[fam.Sample==father_code]['Carrier_Status'])==1:
		father_return = 'FC'
	elif int(fam[fam.Sample==father_code]['Carrier_Status'])==0:
		father_return = 'FNC'
	else:
		father_return = 'F'

	# Check if parents have the same variant
	other_samps = overlap_samps(row['variant_id'])

	if mother_code in other_samps:
		mother = True
	if father_code in other_samps:
		father = True

	# Return inheritance
	if mother and father:
		return 'both'
	elif mother:
		return mother_return
	elif father:
		return father_return

	# If CNV wasn't found in the final output for either person, look into the output files for each CNV caller for each person
	chr = row['Chromosome']
	start = row['Start']
	end = row['End']
	length = row['Length']
	type = row['Type']

	if mother_code != '0' and mother_code in penncnv.PatientID.to_list():
		# Raw PennCNV outputs
		mother = penncnv_lookup(mother_code, chr, start, end, length, type)

	if father_code != '0' and father_code in penncnv.PatientID.to_list():
		# Raw PennCNV outputs
		father = penncnv_lookup(father_code, chr, start, end, length, type)

	if mother and father:
		return 'both'
	elif mother:
		return mother_return
	elif father:
		return father_return
	elif mother_code=='0' or father_code=='0':
		return '.'
	elif mother_code not in penncnv.PatientID.to_list() or father_code not in penncnv.PatientID.to_list():
		return '.'
	else:
		return 'de novo'

calls['inheritance'] = calls.apply(inheritance, axis = 1)

calls.to_csv('call_tables/5_inheritance.txt', sep = '\t', index = False)
