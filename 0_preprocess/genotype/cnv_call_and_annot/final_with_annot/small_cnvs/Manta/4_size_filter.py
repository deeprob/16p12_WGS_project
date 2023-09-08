import pandas as pd


# Iterate through samples and save CNVs to files
samples = pd.read_csv('manta_files.csv')

def split_variants(sample):
	print('Starting '+sample)
	# Read in merged variants
	df = pd.read_csv('bed_files/3_adjacent_filter/'+str(sample)+'_manta_final.bed', sep = '\t')

	# Get the length
	df['Size'] = df.End - df.Pos

	# Add column for sample name - will be needed later
	df['Sample'] = sample

	# Remove CNVs >50kb
	df = df[df.Size < 50000]
	# ALso remove CNVs < 100 bp
	df = df[df.Size >= 100]

	df.to_csv('bed_files/4_size_filter/'+str(sample)+'_size_filter.bed', sep = '\t', index = False)

for sample in samples.Sample.to_list():
	split_variants(sample)
