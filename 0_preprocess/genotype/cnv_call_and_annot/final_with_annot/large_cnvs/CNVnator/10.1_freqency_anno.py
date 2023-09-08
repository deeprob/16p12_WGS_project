import pandas as pd

# Remove samples with an inracohort frequency > 10
# Annotate with microarray frequency

# DELS
df = pd.read_csv('bed_files/9.3_nejm_dels.bed', sep = '\t',
		names = ['Chr', 'Start', 'End', 'Type', 'Name', 'Length', 'Sample', 'NEJM_Name', 'Microarray_count', 'Intracohort_count'])

# Microarray annotations
# 18572 unique individuals in file
# Number taken from Matt's script on Dropbox: Dropbox\16p12.2 project\Human patients project\WGS paper\6_Microarray calls\Analysis\Final calls\control_cnv.py
df['microarray_freq'] = df.Microarray_count/18572
df = df[df.Intracohort_count <= 10]
df.to_csv('bed_files/10_nejm_dels.bed', sep = '\t', index = False)

# DUPS
df = pd.read_csv('bed_files/9.3_nejm_dups.bed', sep = '\t',
		names = ['Chr', 'Start', 'End', 'Type', 'Name', 'Length', 'Sample', 'NEJM_Name', 'Microarray_count', 'Intracohort_count'])

df['microarray_freq'] = df.Microarray_count/18572
df = df[df.Intracohort_count <= 10]
df.to_csv('bed_files/10_nejm_dups.bed', sep = '\t', index = False)
