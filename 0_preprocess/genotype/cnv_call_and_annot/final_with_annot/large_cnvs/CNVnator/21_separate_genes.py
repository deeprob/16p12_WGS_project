import pandas as pd

# Separate calls into gene-level calls
calls = pd.read_csv('final_calls/final_calls.txt', sep = '\t')

# Gene annotations - need these for later
genes = pd.read_csv('/data5/16p12_WGS/structural_variants/sv_caller_postprocessing/2022_02_14/gencode.v19.parsed.genes.csv')
print(genes.columns)
print(genes)

gene_df = pd.DataFrame(columns = ['Sample', 'Gene_ID', 'Gene_Symbol', 'Type', 'call_ID', 'inheritance'])
for idx, row in calls.iterrows():
	# Get call information
	sample = row['Sample']
	gene_ids = row['gene_ids'].split(' ')
	gene_symbols = row['gene_names'].split(' ')
	type = row['Type']
	ID = row['variant_id']
	inheritance = row['inheritance']
    nejm = row['NEJM_Name']

	out_lines = []

	# NOTE: There are some cases where a gene symbol has multiple IDs!!
	# Remove version number from gene IDs
	gene_ids2 = [i.split('.')[0] for i in gene_ids]

	# We only want to report the unique gene symbols
	id_dict = {}
	for n, symbol in enumerate(gene_symbols):
		if symbol not in id_dict:
			id_dict[symbol] = [gene_ids2[n]]
		else:
			id_dict[symbol].append(gene_ids2[n])

	for key in id_dict.keys():
		out_lines.append([sample, ' '.join(id_dict[key]), key, type, nejm, ID, inheritance])

	# Add to dataframe
	gene_df = pd.concat([gene_df, pd.DataFrame(out_lines, columns = ['Sample', 'Gene_ID', 'Gene_Symbol', 'Type', 'NEJM', 'call_ID', 'inheritance'])], axis = 0, ignore_index = True)

print(gene_df)

gene_df.to_csv('final_calls/calls_by_gene.txt', sep = '\t', index = False)
