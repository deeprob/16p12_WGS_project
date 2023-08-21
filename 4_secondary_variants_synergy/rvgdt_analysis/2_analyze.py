#!/usr/bin/env python3


# changing Jianyu's analysis from R to python

DESCRIPTION = "Analyze RVGDT results to find the consensus significant genes across several runs"


import os
import pandas as pd
import seaborn as sns

sns.set(font_scale = 1.75)

def read_rvgdt_output(filename, p, r):
	gene_table = pd.read_csv(filename, sep="\t", header=None, usecols=[0, 2], 
		names=["gene", f"p-value-{r}"], dtype={"gene": str, "p-value": float}
	)
	gene_table["gene"] = gene_table["gene"].str.replace(f".{p}", "", regex=False)
	gene_table = gene_table.set_index("gene")
	return gene_table


phenotypes = ['Child_behav', 'Child_congenital','Child_craniofacial','Child_ID_DD','Child_nervous_system','Child_psych']
runs = list(map(str, range(1, 6)))

rvgdt_results_dir = "/data5/deepro/wgs_16p/rvgdt/data/results"
gene_pval_series = dict()
sig_gene_names = set()

for p in phenotypes:
	filenames = [os.path.join(rvgdt_results_dir, r, f"rvgdt_output_{p}.txt") for r in runs]
	gene_tables = [read_rvgdt_output(fn, p, r) for fn,r in zip(filenames, runs)]
	gene_meta_table = pd.concat(gene_tables, axis=1)
	gene_meta_table_sig = gene_meta_table.loc[(gene_meta_table<0.05).all(axis=1)]
	# append sig gene names to set
	sig_gene_names.update(list(gene_meta_table_sig.index))
	# take the median values of all genes and append in list
	gene_pval_series[p] = gene_meta_table.median(axis=1)
	# save the meta sig table
	meta_filename = os.path.join(rvgdt_results_dir, "tables", f"meta_{p}.csv")
	os.makedirs(os.path.dirname(meta_filename), exist_ok=True)
	gene_meta_table_sig.to_csv(meta_filename)
	

# save a table that contains the gene name and significant value of genes 
# which were significant across all five runs in at least one phenotype
gene_pval_df = pd.DataFrame(gene_pval_series)
gene_pval_df_sig = gene_pval_df.loc[sorted(list(sig_gene_names)), :]
meta_sig_filename = os.path.join(rvgdt_results_dir, "tables", f"meta_sig.csv")
os.makedirs(os.path.dirname(meta_sig_filename), exist_ok=True)
gene_pval_df_sig.to_csv(meta_sig_filename)

# create a clustermap in seaborn
g = sns.clustermap(
	gene_pval_df_sig, 
	row_cluster=False, linewidths=0.25, figsize=(20,30), yticklabels=True,
	cmap=sns.dark_palette("#69d", reverse=True, as_cmap=True), 
	cbar_pos=(1e-10, .2, .02, .4), dendrogram_ratio=0.05)

fig_file = os.path.join(rvgdt_results_dir, "figures", f"meta_sig.pdf")
os.makedirs(os.path.dirname(fig_file), exist_ok=True)
g.savefig(fig_file)
