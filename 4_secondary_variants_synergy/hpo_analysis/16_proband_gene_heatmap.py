#!/bin/python
import pandas as pd
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

matplotlib.rcParams['pdf.fonttype'] = 42

# Make heatmaps showing probands and the genes they have
# These will be later annotated to show phenotypes

# Load in data
# Get probands
cohort_info = pd.read_csv('../../3_cohort information/16p12_All_Participants_v9.csv')
cohort_info = cohort_info[(cohort_info.WGS=='X') & (cohort_info.No_consent_forms!='X')]
pros = cohort_info[cohort_info.Relationship=='P']['Sample'].to_list()

# Phenotypes
pheno = pd.read_csv('../1_HPO_annotations/Annotated_files/phenotype_hpo.csv')[['Sample', 'HPO_ID']]
pheno = pheno[pheno.Sample.isin(pros)]
# Variants
snvs = pd.read_csv('../1_HPO_annotations/Annotated_files/snvs_hpo_anno.csv')[['Sample', 'Gene_symbol', 'Mut_type', 'Gene_id_']]
snvs = snvs[snvs.Sample.isin(pros)]
cnvs = pd.read_csv('../1_HPO_annotations/Annotated_files/cnvs_hpo_anno.csv')[['Sample', 'Gene_Symbol',  'Type', 'Gene_ID']]
cnvs = cnvs[cnvs.Sample.isin(pros)]
# Gene-Phenotype map
gp_map = pd.read_csv('../1_HPO_annotations/Analysis_files/2_ensembl2hpo.csv')

# Get all genes with variants in probands
pro_genes = snvs.Gene_id_.to_list() + cnvs.Gene_ID.to_list()
pro_genes = list(set(pro_genes))
pro_genes.sort()

# Get a list of all probands with variant data
var_pros = [i for i in pros if i in snvs.Sample.to_list() and i in cnvs.Sample.to_list() and i in pheno.Sample.to_list()]

# Make a starting dataframe with phenotypes as rows and probands as columns
gene_df = pd.DataFrame(index = pro_genes, columns = var_pros)

# Make heatmaps
pdf = PdfPages('Figures/2_proband_gene.pdf')
##############################################
# Proportion of phenotypes explained by gene #
##############################################

def prop_explained(col):
    # Get sample ID
    pro = col.name
    
    # Get genes proband has
    genes = []
    dfs = [snvs, cnvs]
    gene_cols = ['Gene_id_', 'Gene_ID']
    for i in range(2):
        df = dfs[i]
        gene_col = gene_cols[i]
        genes += df[df.Sample==pro][gene_col].to_list()
    genes = list(set(genes))
    
    # Get proband phenotypes
    pheno_lst = pheno[pheno.Sample==pro]['HPO_ID'].to_list()[0]
    phenos = [i.strip() for i in pheno_lst.split(';')]
    
    output = []
    for row in pro_genes:
        # Check if proband has gene
        if row not in genes:
            output.append(np.nan)
            continue
        # Check if gene has any HPO terms
        if row not in gp_map.ENSEMBL_ID.to_list():
            output.append(0)
            continue
        # Get HPO terms for gene
        gene_hpo = gp_map[gp_map.ENSEMBL_ID==row]['HPO_IDs'].to_list()[0].split(';')
        
        # Get proportion of proband phenotypes explained by gene HPO
        overlap = [i for i in phenos if i in gene_hpo]
        
        output.append(float(len(overlap))/len(phenos))
    return pd.Series(output, index = col.index)
prop_df = gene_df.apply(lambda col: prop_explained(col), axis = 0)

# Change the row names to gene symbols
snv_symbol = snvs[['Gene_symbol', 'Gene_id_']]
cnv_symbol = cnvs[['Gene_Symbol', 'Gene_ID']]
cnv_symbol.columns = ['Gene_symbol', 'Gene_id_']
symbol_df = pd.concat([snv_symbol, cnv_symbol])
del snv_symbol
del cnv_symbol
symbol_df.set_index(symbol_df.Gene_id_, inplace = True)
symbol_dict = symbol_df.Gene_symbol.to_dict()
prop_df.index = prop_df.index.map(symbol_dict)
# Save to file
out_df = prop_df.fillna('.')
out_df.to_csv('Result_tables/2_proband_gene_proportion.csv')

# Plot heatmap

# Function for making square map
def fixedWidthClusterMap(dataFrame, colors, cellSizePixels=50, cbar = (0, .2, .03, .4), clust=True):
    # Calulate the figure size, this gets us close, but not quite to the right place
    dpi = matplotlib.rcParams['figure.dpi']
    marginWidth = matplotlib.rcParams['figure.subplot.right']-matplotlib.rcParams['figure.subplot.left']
    marginHeight = matplotlib.rcParams['figure.subplot.top']-matplotlib.rcParams['figure.subplot.bottom']
    Ny,Nx = dataFrame.shape
    figWidth = (Nx*cellSizePixels/dpi)/0.8/marginWidth
    figHeigh = (Ny*cellSizePixels/dpi)/0.8/marginHeight

    # do the actual plot
    grid = sns.clustermap(dataFrame, figsize=(figWidth, figHeigh), cmap = colors, cbar_pos = cbar, linewidths = 0.3, linecolor = 'lightgrey',
                            dendrogram_ratio = (0.15, 0.00001), xticklabels = True, yticklabels = True, row_cluster = clust, col_cluster = clust)

    # calculate the size of the heatmap axes
    axWidth = (Nx*cellSizePixels)/(figWidth*dpi)
    axHeight = (Ny*cellSizePixels)/(figHeigh*dpi)

    # resize heatmap
    ax_heatmap_orig_pos = grid.ax_heatmap.get_position()
    grid.ax_heatmap.set_position([ax_heatmap_orig_pos.x0, ax_heatmap_orig_pos.y0, 
                                  axWidth, axHeight])

    # resize dendrograms to match
    ax_row_orig_pos = grid.ax_row_dendrogram.get_position()
    grid.ax_row_dendrogram.set_position([ax_row_orig_pos.x0, ax_row_orig_pos.y0, 
                                         ax_row_orig_pos.width, axHeight])
    ax_col_orig_pos = grid.ax_col_dendrogram.get_position()
    grid.ax_col_dendrogram.set_position([ax_col_orig_pos.x0, ax_heatmap_orig_pos.y0+axHeight,
                                         axWidth, ax_col_orig_pos.height])
    return grid # return ClusterGrid object

# Plot is too large with all genes/probands
# Filter genes for those that explain phenotypes in at least 1 proband
expl_df = prop_df[prop_df > 0].dropna(thresh = 1, axis = 0)
prop_df = prop_df.loc[expl_df.index, expl_df.columns]
prop_df.fillna(-0.1, inplace = True)

# Make plot
colors = ['white', 'black'] + sns.color_palette("mako", 10)[0:9]
cg = fixedWidthClusterMap(prop_df, colors, cellSizePixels = 10)
cg.ax_row_dendrogram.set_visible(False)
cg.ax_col_dendrogram.set_visible(False)
cg.ax_heatmap.set_xticklabels(cg.ax_heatmap.get_xmajorticklabels(), fontsize = 4)
cg.ax_heatmap.set_yticklabels(cg.ax_heatmap.get_ymajorticklabels(), fontsize = 6)
cg.fig.suptitle('Proportions: All Probands, Genes explain >= 1 proband')
pdf.savefig()
plt.close()

# Filter genes for those that explain phenotypes in at least 3 probands
prop_df.replace(-0.1, np.nan, inplace = True)
expl_df = prop_df[prop_df > 0].dropna(thresh = 3, axis = 0)
# Also remove probands who do not have any of these genes
expl_df = expl_df[expl_df > 0].dropna(thresh = 1, axis = 1)
prop_df = prop_df.loc[expl_df.index, expl_df.columns]
prop_df.fillna(-0.1, inplace = True)

# Make plot
cg2 = fixedWidthClusterMap(prop_df, colors, cellSizePixels = 10)
cg2.ax_row_dendrogram.set_visible(False)
cg2.ax_col_dendrogram.set_visible(False)
cg2.ax_heatmap.set_xticklabels(cg2.ax_heatmap.get_xmajorticklabels(), fontsize = 4)
cg2.ax_heatmap.set_yticklabels(cg2.ax_heatmap.get_ymajorticklabels(), fontsize = 6)
cg2.fig.suptitle('Proportions: All Probands, Genes explain >= 3 probands')
pdf.savefig()
plt.close()


###########################
# Type of variant in gene #
###########################

# Codes:
# 1. LOF variant
# 2. Missense variant
# 3. Splice variant
# 4. Deletion
# 5. Duplication
# 6. Multiple variants
def type_explain(col):
    # Get sample ID
    pro = col.name
    
    # Get genes proband has
    genes = []
    dfs = [snvs, cnvs]
    gene_cols = ['Gene_id_', 'Gene_ID']
    for i in range(2):
        df = dfs[i]
        gene_col = gene_cols[i]
        genes += df[df.Sample==pro][gene_col].to_list()
    genes = list(set(genes))
    
    output = []
    for row in pro_genes:
        # Check if proband has gene
        if row not in genes:
            output.append(np.nan)
            continue
        # Get type of variant proband has
        dfs = [snvs, cnvs]
        types = []
        for i in range(2):
            df = dfs[i]
            gene_col = gene_cols[i]
            sub_df = df[(df.Sample==pro) & (df[gene_col]==row)]
            if sub_df.shape[0]==0:
                # Gene has no variant in that dataframe for that proband
                continue
            if i==0:
                for type in ['missense', 'lof', 'splice']:
                    sub_df2 = sub_df[sub_df.Mut_type==type]
                    if sub_df2.shape[0]==0:
                        continue
                    types.append(type)
            elif i==1:
                for type in ['DEL', 'DUP']:
                    sub_df2 = sub_df[sub_df.Type==type]
                    if sub_df2.shape[0]==0:
                        continue
                    types.append(type)
        if len(types) > 1:
            val = 6
        elif types==['lof']:
            val = 1
        elif types==['missense']:
            val = 2
        elif types==['splice']:
            val = 3
        elif types==['DEL']:
            val = 4
        elif types==['DUP']:
            val = 5
        output.append(val)
    return pd.Series(output, index = col.index)
type_df = gene_df.apply(lambda col: type_explain(col), axis = 0)

# Add gene names
type_df.index = type_df.index.map(symbol_dict)
# Organize rows and columns to match last heatmap
row_names = [i.get_text() for i in cg2.ax_heatmap.get_ymajorticklabels()]
col_names = [i.get_text() for i in cg2.ax_heatmap.get_xmajorticklabels()]
type_df = type_df.loc[row_names, col_names]
# Change NA to -100
type_df.fillna(-100, inplace = True)

# Make plots
colors = ['white']*101 + sns.color_palette("tab10")[0:5]
cg3 = fixedWidthClusterMap(type_df, colors, cellSizePixels = 10, clust = False)
cg3.ax_row_dendrogram.set_visible(False)
cg3.ax_col_dendrogram.set_visible(False)
cg3.ax_heatmap.set_xticklabels(cg3.ax_heatmap.get_xmajorticklabels(), fontsize = 4)
cg3.ax_heatmap.set_yticklabels(cg3.ax_heatmap.get_ymajorticklabels(), fontsize = 6)
cg3.fig.suptitle('Types: All Probands, Genes explain >= 3 probands')
pdf.savefig()
plt.close()


# Filter for genes with multiple types of variants
type_df=type_df[type_df.apply(lambda row: len(list(set(row.to_list())))>2, axis=1)]
keep_cols=type_df.columns[type_df.apply(lambda col: len(list(set(col.to_list())))>1, axis=0)]
type_df=type_df[keep_cols]
# Make plots
colors = ['white']*101 + sns.color_palette("tab10")[0:5]
cg3 = fixedWidthClusterMap(type_df, colors, cellSizePixels = 10, clust = False)
cg3.ax_row_dendrogram.set_visible(False)
cg3.ax_col_dendrogram.set_visible(False)
cg3.ax_heatmap.set_xticklabels(cg3.ax_heatmap.get_xmajorticklabels(), fontsize = 4)
cg3.ax_heatmap.set_yticklabels(cg3.ax_heatmap.get_ymajorticklabels(), fontsize = 6)
cg3.fig.suptitle('Types: All Probands, Genes explain >= 3 probands, Genes have multiple variant types')
pdf.savefig()
plt.close()

# Proportion
prop_df=prop_df.loc[type_df.index, type_df.columns]
cg3 = fixedWidthClusterMap(prop_df, ['white', 'black'] + sns.color_palette("mako", 10)[0:9], cellSizePixels = 10, clust = False)
cg3.ax_row_dendrogram.set_visible(False)
cg3.ax_col_dendrogram.set_visible(False)
cg3.ax_heatmap.set_xticklabels(cg3.ax_heatmap.get_xmajorticklabels(), fontsize = 4)
cg3.ax_heatmap.set_yticklabels(cg3.ax_heatmap.get_ymajorticklabels(), fontsize = 6)
cg3.fig.suptitle('Types: All Probands, Genes explain >= 3 probands, Genes have multiple variant types')
pdf.savefig()
plt.close()

pdf.close()
