import pandas as pd

# Remove samples with
# A microarray control frequency > 0.1
# or, an inracohort frequency > 10

# DELS
df = pd.read_csv('bed_files/9.3_dels.bed', sep = '\t',
		names = ['Chr', 'Start', 'End', 'Type', 'Name', 'Length', 'Sample', 'Microarray_count', 'Intracohort_count'])

# Microarray filter
# 18572 unique individuals in file
# Number taken from Matt's script on Dropbox: Dropbox\16p12.2 project\Human patients project\WGS paper\6_Microarray calls\Analysis\Final calls\control_cnv.py
df['microarray_freq'] = df.Microarray_count/18572
df = df[df.microarray_freq <= 0.001]

# Intracohort filter
df = df[df.Intracohort_count <= 10]

df.to_csv('bed_files/10_dels_frequency_filter.bed', sep = '\t', index = False)

# DUPS
df = pd.read_csv('bed_files/9.3_dups.bed', sep = '\t',
		names = ['Chr', 'Start', 'End', 'Type', 'Name', 'Length', 'Sample', 'Microarray_count', 'Intracohort_count'])

df['microarray_freq'] = df.Microarray_count/18572
df = df[df.microarray_freq <= 0.001]
df = df[df.Intracohort_count <= 10]

df.to_csv('bed_files/10_dups_frequency_filter.bed', sep = '\t', index = False)
