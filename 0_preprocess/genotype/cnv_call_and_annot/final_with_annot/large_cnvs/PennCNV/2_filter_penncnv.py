import pandas as pd

# Filter PennCNV calls based on frequency and gene number
calls = pd.read_csv('../../2022_02_14/PennCNV_calls_combined_final.csv')
print(calls)

# Merge with NEJM CNVs
nejm = pd.read_csv('bed_files/1_nejm_filter.bed', sep = '\t', names = ['Chromosome', 'Start', 'End', 'Type', 'PatientID', 'NEJM_Chr', 'NEJM_Start', 'NEJM_End', 'NEJM_Name'])
print(nejm)

calls = pd.merge(calls, nejm[['Chromosome', 'Start', 'End', 'Type', 'PatientID','NEJM_Name']], on = ['Chromosome', 'Start', 'End', 'Type', 'PatientID'], how = 'left')
print(calls)

# Replace empty NEJM names with .
calls.loc[calls.NEJM_Name.isnull(), "NEJM_Name"] = '.'
print(calls)
print(calls.NEJM_Name.value_counts())

# Intracohort frequency filter
calls = calls[calls.intracohort_count <= 10]
# gnomAD frequency filter (do not apply this to NEJM CNVs)
calls = calls[(calls.gnomad_freq_overall <= 0.001) | (calls.NEJM_Name != '.')]
print(calls.NEJM_Name.value_counts())

# Filter for SegDups and Centromere/Telomere
# 50% overlap - use the present annotations
# Do not apply this to NEJM CNVs
calls = calls[(calls["%SD"] <= 0.5) | (calls.NEJM_Name != '.')]
calls = calls[(calls["%CenTel"] <= 0.5) | (calls.NEJM_Name != '.')]

# Also do some reformatting
calls.Type = calls.Type.apply(lambda row: row.upper())

# Add a variant ID
calls['variant_id'] = calls.Chromosome+'_'+calls.Start.astype(str)+'_'+calls.End.astype(str)+'_'+calls.Type+'_'+calls.PatientID

# Save to file
calls.to_csv('call_tables/2_filter.txt', sep = '\t', index = False)
