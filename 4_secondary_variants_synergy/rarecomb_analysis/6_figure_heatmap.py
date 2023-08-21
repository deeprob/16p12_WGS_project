#!/bin/python3


import pandas as pd
import sys
import numpy as np
import subprocess

# libraries related to plotting
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
sns.set_style({'font.family':'sans-serif', 'font.sans-serif':'Arial'})
from matplotlib.backends.backend_pdf import PdfPages


subprocess.run('mkdir figures', shell=True)


df2 = pd.read_csv(f'statistics/5_combine/combolen2.csv')
df3 = pd.read_csv(f'statistics/5_combine/combolen3.csv')


genes = df2['Item_1'].to_list()
genes = genes + df2['Item_2'].to_list()
genes = genes + df3['Item_1'].to_list()
genes = genes + df3['Item_2'].to_list()
genes = genes + df3['Item_3'].to_list()


genes = pd.Series(genes)
gene_counts = genes.value_counts()
gene_counts = pd.DataFrame(gene_counts)
gene_counts['Gene'] = gene_counts.index.to_series()
gene_counts['Count'] = gene_counts[0]
genes = gene_counts['Gene'].to_list()


# dfplot = pd.DataFrame(index=genes)
dfplot = pd.DataFrame(columns=['Phenotype'] + genes)
current_row = 0
for i, row in df2.iterrows():
	gene1 = row['Item_1']
	gene2 = row['Item_2']
	phenotype = row['Phenotype']
	fdr = row['Case_Adj_Pval_bonf']
	dfplot.loc[current_row] = 0
	dfplot.loc[current_row, 'Phenotype'] = f'{phenotype} ({fdr})'
	dfplot.loc[current_row, gene1] = 1
	dfplot.loc[current_row, gene2] = 1
	current_row = current_row + 1


for i, row in df3.iterrows():
	gene1 = row['Item_1']
	gene2 = row['Item_2']
	gene3 = row['Item_3']
	phenotype = row['Phenotype']
	fdr = row['Case_Adj_Pval_bonf']
	dfplot.loc[current_row] = 0
	dfplot.loc[current_row, 'Phenotype'] = f'{phenotype} ({fdr})'
	dfplot.loc[current_row, gene1] = 1
	dfplot.loc[current_row, gene2] = 1
	dfplot.loc[current_row, gene3] = 1
	current_row = current_row + 1


dfplot = dfplot.sort_values('Phenotype')
dfplot = dfplot.set_index('Phenotype', drop=True)



for col in dfplot.columns:
	dfplot[col] = dfplot[col].astype(float)




pdf = PdfPages('figures/heatmap_rarecomb.pdf')

fig = plt.figure(figsize=(15,5))

g = sns.heatmap(data=dfplot, fmt='', square=True, cmap='bwr', center=0, linecolor='black', linewidths=0.1, cbar=False)
pdf.savefig(fig, bbox_inches='tight')


pdf.close()



