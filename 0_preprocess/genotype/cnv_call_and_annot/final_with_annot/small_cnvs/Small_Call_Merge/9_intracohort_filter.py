import pandas as pd

# Remove samples with an inracohort frequency > 10

# DELS
df = pd.read_csv('bed_files/8_dels.bed', sep = '\t',
		names = ['Chr', 'Start', 'End', 'Type', 'Name', 'Length', 'Sample', 'Intracohort_count'])
# Intracohort filter
df = df[df.Intracohort_count <= 10]
df.to_csv('bed_files/9_dels_intracohort_filter.bed', sep = '\t', index = False)

# DUPS
df = pd.read_csv('bed_files/8_dups.bed', sep = '\t',
		names = ['Chr', 'Start', 'End', 'Type', 'Name', 'Length', 'Sample', 'Intracohort_count'])

df = df[df.Intracohort_count <= 10]
df.to_csv('bed_files/9_dups_intracohort_filter.bed', sep = '\t', index = False)
