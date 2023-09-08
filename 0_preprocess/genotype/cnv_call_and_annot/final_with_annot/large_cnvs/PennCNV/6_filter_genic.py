import pandas as pd

# Filter calls for only those that overlap genes
calls = pd.read_csv('call_tables/5_inheritance.txt', sep = '\t')
print(calls)

# Remove those that overlap no genes
calls = calls[calls.gene_ids!='.']

# Remove unneeded columns
calls = calls[['PatientID', 'Chromosome', 'Start', 'End', 'Type', 'Length', 'NumSNPs', 'StartSNP', 'EndSNP', '%SD', '%CenTel', 'intracohort_count', 'gnomad_freq_overall',
		'NEJM_Name', 'variant_id', 'gene_ids', 'gene_names', 'inheritance']]

# Rename some columns
calls.columns = ['Sample', 'Chr', 'Start', 'End', 'Type', 'Length', 'NumSNPs', 'StartSNP', 'EndSNP', '%SD', '%CenTel', 'intracohort_count', 'gnomad_freq_overall',
		'NEJM_Name', 'variant_id', 'gene_ids', 'gene_names', 'inheritance']

print(calls)

# Save to file
calls.to_csv('call_tables/6_genic_filter.txt', sep = '\t', index = False)
