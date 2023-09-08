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
calls = pd.read_csv('bed_files/14_gene_filter.bed', sep = '\t')
lookup = pd.read_csv('bed_files/15_inheritance_lookup.bed', sep = '\t',
	names = ['Chr', 'Start', 'End', 'Type', 'Name', 'Length', 'Sample', 'Intracohort_count', 'gnomADSV_AF', 'gnomADSV_AFfilter', 'gene_ids', 'gene_names', 'variant_id',
		'Chr2', 'Start2', 'End2', 'Type2', 'Name2', 'Length2', 'Sample2', 'Intracohort_count2', 'gnomADSV_AF2', 'gnomADSV_AFfilter2', 'gene_ids2', 'gene_names2', 'variant_id2'])

fam = pd.read_csv('../../2022_02_14/Feb_2022.fam', sep = '\t', names = ['Family', 'Sample', 'Father', 'Mother', 'Sex', 'Carrier_Status'])

cnv_samps = pd.read_csv('small_cnv_samples.list', names = ['Sample'])

# Caller output file locations
cnvnator = pd.read_csv('../../2022_02_14/cnvnator/cnvnator_files.csv')
manta =  pd.read_csv('../../2022_02_14/manta/manta_files.csv')
delly = pd.read_csv('../../2022_02_14/delly/delly_files.csv')
lumpy = pd.read_csv('../lumpy/lumpy_files.csv')

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

# Use bcftools to search through the merged output from the callers to search for calls that didn't pass QC
def rawfile_lookup(file, chr, start, end, length, svtype):
	# Due to 50% reciprocal overlap requirement, any overlapping CNV can't start more than LENGTH away from the START of the query CNV
	# It also can't start any more than half way through the query CNV
	search_start = start-length
	search_end = start+(length/2)

	command = "bcftools view -r %s:%i-%i -t %s:%i-%i -i 'INFO/SVTYPE=\"%s\"' %s -H" % (chr, search_start, search_end, chr, search_start, search_end, svtype, file)
	match_lines = subprocess.run(command, capture_output = True, shell = True).stdout.decode()

	if len(match_lines) == 0:
		return False

	match_list = [i.split('\t') for i in match_lines.split('\n') if i!='']
	# Convert to df
	match_df = pd.DataFrame(match_list, columns = ['Chr','Pos', 'ID', 'Ref', 'Alt', 'Qual', 'FILTER', 'INFO', 'FORMAT', 'FORMAT_VALUES'])
	# Separate out end and length
	match_df['End'] = match_df.INFO.str.split('END=', expand = True)[1].str.split(';', expand = True)[0]
	match_df['Pos'] = pd.to_numeric(match_df['Pos'])
	match_df['End'] = pd.to_numeric(match_df['End'])

	match_df['Length'] = match_df.End - match_df.Pos

	# Get only calls that meet 50% reciprocal overlap
	match_df = match_df[match_df.apply(lambda row: recip_overlap(start, end, row['Pos'], row['End'], length, row['Length']), axis = 1)]

	if match_df.shape[0] > 0:
		return True

	return False

# Look through merged files from each caller for overlap
def mergefile_lookup(file,chr, start, end, length, svtype):
	# Read in parent file as dataframe
	df = pd.read_csv(file, sep = '\t')
	# Unify column names
	df.columns = ['Chr', 'Pos', 'End', 'Type', 'ID']
	# Unify SV type name
	type_dict = {'<DEL>':'DEL', '<DUP>':'DUP', 'DEL':'DEL', 'DUP':'DUP'}
	df['Type'] = df.Type.map(type_dict)

	# Filter line
	search_start = start-length
	search_end = start+(length/2)
	df = df[(df.Chr==chr) & (df.Type==svtype) & (df.Pos>=search_start) & (df.Pos<=search_end)]

	# Add length
	df['Length'] = df.End - df.Pos

	# Get 50% reciprocal overlap calls
	df = df[df.apply(lambda row: recip_overlap(start, end, row['Pos'], row['End'], length, row['Length']), axis = 1)]

	if df.shape[0] > 0:
		return True

	return False

def inheritance(row, father = False, mother = False):
	iid = row['Sample']

	# Get family information from FAM file
	mother_code = fam[fam.Sample==iid]['Mother'].to_string(index = False, header = False).strip()
	father_code = fam[fam.Sample==iid]['Father'].to_string(index = False, header = False).strip()

	# Check if we have CNV calls for the parents
	if mother_code == '0' and father_code == '0':
		return '.'
	if mother_code not in cnv_samps.Sample.to_list() and father_code not in cnv_samps.Sample.to_list():
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
	chr = row['Chr']
	start = row['Start']
	end = row['End']
	length = row['Length']
	type = row['Type']

	caller_list_files = [cnvnator, manta, delly, lumpy]
	caller_names = ['cnvnator', 'manta', 'delly', 'lumpy']
	if mother_code != '0':
		for i, caller in enumerate(caller_list_files):
			if mother_code in caller.Sample.to_list():
				name = caller_names[i]
				# Filtered outputs
				mother_file = "compressed_caller_outputs/"+name+"/"+mother_code+"."+name+".filtered.vcf.gz"
				mother = rawfile_lookup(mother_file, chr, start, end, length, type)

				if mother:
					break

				# Merged outputs
				if name=='lumpy':
					mother_merge = "../lumpy/bed_files/3_adjacent_filter/"+mother_code+"_lumpy_final.bed"
				elif name=='cnvnator':
					mother_merge = "../../2022_02_14/cnvnator/tables/2_adjacent_filter/"+mother_code+"_cnvnator_final.bed"
				else:
					mother_merge = "../../2022_02_14/"+name+"/bed_files/3_adjacent_filter/"+mother_code+"_"+name+"_final.bed"
				mother = mergefile_lookup(mother_merge, chr, start, end, length, type)

				if mother:
					break


	if father_code != '0':
		for i, caller in enumerate(caller_list_files):
			if father_code in caller.Sample.to_list():
				name = caller_names[i]
				# Filtered outputs
				father_file = "compressed_caller_outputs/"+name+"/"+father_code+"."+name+".filtered.vcf.gz"
				father = rawfile_lookup(father_file, chr, start, end, length, type)

				if father:
					break

				# Merged outputs
				if name=='lumpy':
					father_merge = "../lumpy/bed_files/3_adjacent_filter/"+father_code+"_lumpy_final.bed"
				elif name=='cnvnator':
					father_merge = "../../2022_02_14/cnvnator/tables/2_adjacent_filter/"+father_code+"_cnvnator_final.bed"
				else:
					father_merge = "../../2022_02_14/"+name+"/bed_files/3_adjacent_filter/"+father_code+"_"+name+"_final.bed"
				father = mergefile_lookup(father_merge, chr, start, end, length, type)

				if father:
					break

	if mother and father:
		return 'both'
	elif mother:
		return mother_return
	elif father:
		return father_return
	elif mother_code=='0' or father_code=='0':
		return '.'
	elif mother_code not in cnv_samps.Sample.to_list() or father_code not in cnv_samps.Sample.to_list():
		return '.'
	else:
		return 'de novo'

calls['inheritance'] = calls.apply(inheritance, axis = 1)

calls.to_csv('bed_files/16_inheritance.bed', sep = '\t', index = False)
