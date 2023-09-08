import pandas as pd

# Limit calls to only those we are using for WGS analysis
calls = pd.read_csv('call_tables/8_update_inheritance.txt', sep = '\t')
print(calls.Sample.value_counts())

table = pd.read_csv('/data5/16p12_WGS/structural_variants/sv_caller_postprocessing/2022_02_14/16p12_All_Participants_v5.csv')
valid_samp = table[table.WGS=='X']['Sample'].to_list()

calls = calls[calls.Sample.isin(valid_samp)]
print(calls.Sample.value_counts())

# Also make sure all samples passed microarray
valid_micro = table[table.Microarray.isin(['2019', '2017'])]['Sample'].to_list()
calls = calls[calls.Sample.isin(valid_micro)]
print(calls.Sample.value_counts())

calls.to_csv('call_tables/9_finalize_calls.txt', sep = '\t', index = False)

