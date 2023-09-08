import pandas as pd

# Filter the calls using the gnomAD SV AF
# Dels
lookup = pd.read_csv('bed_files/11.1_nejm_dels.bed', sep = '\t')
dels = pd.read_csv('bed_files/10_nejm_dels.bed', sep = '\t')

# There are some CNVs that have 50% reciprocal overlap with more than one SV in gnomAD
# I will filter based on the highest AF, but report both
def get_AFs(s):
	chr = s.Chr
	start = s.Start
	end = s.End
	# Find lines in the lookup table with the same location
	overlap_lines = lookup[(lookup.Chr==chr) & (lookup.Start==start) & (lookup.End==end)]
	# Get allele frequencies
	afs = list(overlap_lines.AF.unique())
	afs.sort()
	s['gnomADSV_AF'] = ' '.join(afs)

	if afs==['.']:
		s['gnomADSV_AFfilter'] = True
	elif max([float(i) for i in afs if i!='.']) <= 0.001:
		s['gnomADSV_AFfilter'] = True
	else:
		s['gnomADSV_AFfilter'] = False

	return s

dels = dels.apply(get_AFs, axis = 1)
dels.to_csv('bed_files/12.1_nejm_dels.bed', index = False, sep = '\t')

# Dups
lookup = pd.read_csv('bed_files/11.1_nejm_dups.bed', sep = '\t')
dups = pd.read_csv('bed_files/10_nejm_dups.bed', sep = '\t')
dups = dups.apply(get_AFs, axis = 1)
dups.to_csv('bed_files/12.1_nejm_dups.bed', index = False, sep = '\t')
