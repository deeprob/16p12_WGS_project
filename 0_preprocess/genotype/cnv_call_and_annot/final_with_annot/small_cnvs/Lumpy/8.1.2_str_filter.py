import pandas as pd
import sys

# Read in CNVnator calls
cnvs = pd.read_csv(sys.argv[1], sep = '\t')
print(cnvs)

sample = sys.argv[1].split('/')[2].split('_')[0]

# Read in STR regions
strs = pd.read_csv('../../2022_02_14/hg19_ver13_1.bed', sep = '\t', names = ['Chr', 'Start', 'End', 'motif_length', 'motif'])
print(strs)

# Remove any CNV with a breakpoint in an STR region
def str_remove(chr, start, end):
	str_chr = strs[strs.Chr == chr]
	# Get rows in STR file that overlap a breakpoint
	str_breakpoint = str_chr[((str_chr.Start <= start) & (str_chr.End>=start)) | ((str_chr.Start <= end) & (str_chr.End >= end))]
	# If any STRs overlap breakpoints, remove it from output
	if str_breakpoint.shape[0] == 0:
		return True
	else:
		return False

cnvs = cnvs[cnvs.apply(lambda row: str_remove(row['Chr'], row['Pos'], row['End']), axis = 1)]

# Write to file
cnvs.to_csv('bed_files/8_str_overlap/'+sample+'_8.3_all_str_filter.bed', sep = '\t', index = False, header = False)
