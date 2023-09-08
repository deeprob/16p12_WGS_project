import pandas as pd

# Split CNVs into Dels and Dups for easy frequency filtering
df = pd.read_csv('bed_files/6_combine_samples.bed', sep = '\t')

del_df = df[df.Type=='DEL']
dup_df = df[df.Type=='DUP']

del_df.to_csv('bed_files/7_dels.bed', sep = '\t', index = False)
dup_df.to_csv('bed_files/7_dups.bed', sep = '\t', index = False)
