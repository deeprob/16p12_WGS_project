#!/bin/python
import pandas as pd

# Reformat the gene-HPO term links
df = pd.read_csv('Analysis_files/genes_to_phenotype.txt', sep = '\t', header = None, skiprows = 1,
    names = ['entrez-gene-id', 'entrez-gene-symbol', 'HPO-Term-ID', 'HPO-Term-Name', 'Frequency-Raw', 'Frequency-HPO', 'Additional Info from G-D source', 'G-D source', 'disease-ID for link'])

# Save output to a new dataframe
out_df = df[['entrez-gene-id', 'entrez-gene-symbol']].copy()
out_df.drop_duplicates(inplace = True)
print(out_df)

# Add HPO annotations
def add_hpo(gene_id):
    rel_rows = df[df['entrez-gene-id']==gene_id]
    hpo = list(rel_rows['HPO-Term-ID'].unique())
    return ';'.join(hpo)

out_df['HPO_IDs'] = out_df['entrez-gene-id'].apply(add_hpo)

# Save to file
out_df.to_csv('Analysis_files/1_hpo_gene_annotations.csv', index = False)