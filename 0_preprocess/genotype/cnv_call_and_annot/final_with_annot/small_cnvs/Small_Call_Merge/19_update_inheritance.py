import pandas as pd

# For de novo small CNVs, we decided to vizualize with samplot and only take calls that look real
# For any call that did not seem de novo, we will instead call as unknown inheritance (".")

# Visualizations were done with sctipt 17_denovo_check.sh and uploaded to Dropbox (Dropbox\16p12.2 project\Human patients project\WGS paper\7_Structural variants\samplots\small_calls\denovo_plots)
# I (Corrine) manually looked through all calls and chose those that looked de novo
# My notes are in the file small_DN_calls.txt (on Dropbox at:Dropbox\16p12.2 project\Human patients project\WGS paper\7_Structural variants\samplots\small_calls\small_DN_calls.txt)

# Merge hand-annotations and scipt-annotations
hand_anno = pd.read_csv('small_DN_calls.txt', sep = '\t')
print(hand_anno)
calls = pd.read_csv('bed_files/16_inheritance.bed', sep = '\t')
print(calls)

calls = pd.merge(calls, hand_anno[['variant_id', 'Actual inheritance']], on = 'variant_id', how = 'left')
print(calls.inheritance.value_counts())

# Update inheritance
def update_inherit(row):
	inherit = row['inheritance']
	inherit2 = row['Actual inheritance']

	# Update "both" to unknown
	if inherit=='both':
		return '.'

	# No updated inheritance
	if inherit2!=inherit2:
		return inherit

	# Cases we need to update from hand annotations
	if inherit=='de novo':
		# If hand-annotation is also de novo
		if inherit2=='de novo':
			return 'de novo'
		# Otherwise return unknown
		else:
			return '.'

calls['inheritance'] = calls.apply(update_inherit, axis = 1)
print(calls.inheritance.value_counts())

# Remove extra column
calls.drop('Actual inheritance', axis = 1, inplace = True)

# Reorder columns and save to file
calls = calls[['Sample', 'Chr', 'Start', 'End', 'Type', 'Name', 'Length', 'Intracohort_count', 'gnomADSV_AF', 'gene_ids', 'gene_names', 'variant_id', 'inheritance']]
calls.to_csv('final_calls/final_calls.txt', sep = '\t', index = False)

