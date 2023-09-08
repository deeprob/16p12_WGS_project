import pandas as pd

# March 8 2022
# After discussion with Santhosh and Matt, we decided to no longer use de novo calls from CNVnator
# We will keep inherited calls, but all de novo will be converted to "unknown"
# We will only use microarray for reporting de novo calls

calls = pd.read_csv('bed_files/17_inheritance.bed', sep = '\t')

def update_inherit(inherit):
	# Convert de novo to unknown
	if inherit == 'de novo':
		return '.'
	# Also convert "both" to unknown
	if inherit == 'both':
		return '.'
	return inherit

calls['inheritance'] = calls.inheritance.apply(update_inherit)

# Also update Type
type_dict = {'<DEL>':'DEL', '<DUP>':'DUP'}
calls['Type'] = calls.Type.map(type_dict)

# And NEJM Name
def update_nejm(nejm):
	# If there is no name
	if nejm!=nejm:
		return '.'
	return nejm

# Filter for only genic calls
calls = calls[calls.gene_ids!='.']

# Reorder columns and save to new file
# Reorder columns and save to file
calls = calls[['Sample', 'Chr', 'Start', 'End', 'Type', 'Name', 'Length', 'Intracohort_count', 'Microarray_count', 'microarray_freq', 'gnomADSV_AF',
		'NEJM_Name', 'gene_ids', 'gene_names', 'variant_id', 'inheritance']]
calls.to_csv('final_calls/final_calls.txt', sep = '\t', index = False)
