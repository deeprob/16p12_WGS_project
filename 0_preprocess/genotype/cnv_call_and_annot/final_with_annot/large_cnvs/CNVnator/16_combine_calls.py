import pandas as pd

# Combined NEJM and other calls into one callset
# Make sure to remove duplicates
large = pd.read_csv('bed_files/15_gene_anno.bed', sep = '\t')
nejm = pd.read_csv('bed_files/15_nejm_gene_anno.bed', sep = '\t')

combined = pd.concat([large, nejm], axis = 0, ignore_index = True)

# Remove duplicates
combined['variant_id'] = combined['Chr']+'_'+combined['Start'].astype(str)+'_'+combined['End'].astype(str)+'_'+combined['Sample']
combined.drop_duplicates(subset='variant_id', keep = 'last', inplace = True)

# While I'm here . . . remove chromosome Y calls
combined = combined[combined.Chr!='chrY']

# Save to new file
combined.to_csv('bed_files/16_all_calls_combined.bed', sep = '\t', index = False)
