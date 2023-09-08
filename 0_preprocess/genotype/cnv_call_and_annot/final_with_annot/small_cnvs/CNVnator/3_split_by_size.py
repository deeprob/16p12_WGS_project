import pandas as pd


# Iterate through samples and save CNVs to files
samples = pd.read_csv('cnvnator_files.csv')

def split_variants(sample):
	print('Starting '+sample)
	# Read in merged variants
	df = pd.read_csv('tables/2_adjacent_filter/'+str(sample)+'_cnvnator_final.bed', sep = '\t')

	# Get the length
	df['Size'] = df.End - df.Pos

	# Add column for sample name - will be needed later
	df['Sample'] = sample

	#Separate >=50kb CNVs and <50kb CNVs
	large_df = df[df.Size >= 50000]
	small_df = df[df.Size < 50000]
	# Also remove CNVs < 100bp
	small_df = small_df[small_df.Size >= 100]

	large_df.to_csv('large_cnv_processing/bed_files/3_split_by_size/'+str(sample)+'_cnvnator_large.bed', sep = '\t', index = False)
	small_df.to_csv('small_cnv_processing/bed_files/3_split_by_size/'+str(sample)+'_cnvnator_small.bed', sep = '\t', index = False)
	print('Large CNVs: '+str(large_df.shape[0]))
	print('Small CNVs: '+str(small_df.shape[0]))

for sample in samples.Sample.to_list():
	split_variants(sample)
