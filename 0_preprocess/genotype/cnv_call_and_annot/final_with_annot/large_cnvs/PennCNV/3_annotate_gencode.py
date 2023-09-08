
import pandas as pd

# Annotate gencode genes using gencode.v19.parsed.exons.csv

# gencode.v19.parsed.exons.csv was originally made by Anastasia
# Anastasia's parsing scripts are here: Dropbox/16p12.2 project/Human patients project/WGS paper/5_SNV pipelines-annotations/GENCODE annotations
gencode = pd.read_csv('../../2022_02_14/gencode.v19.parsed.exons.csv')
gencode = gencode[gencode.gene_type=='protein_coding']

# Large CNVs
cnv = pd.read_csv('call_tables/2_filter.txt', sep = '\t')

# Get genes
def get_genes(s):
	chr = s.Chromosome
	start = s.Start
	end = s.End

	# Get exons from gencode file that are at least partially contained within the CNV length
	gene_df = gencode[(gencode.Chrom==chr) & (((gencode.Start>=start) & (gencode.Start<=end)) | ((gencode.End>=start) & (gencode.End<=end)))]

	if gene_df.shape[0]>0:
		# Get a list of gene IDs and gene names
		gene_info = gene_df[['gene_id', 'gene_name']].drop_duplicates()
		s['gene_ids'] = ' '.join(gene_info.gene_id)
		s['gene_names'] = ' '.join(gene_info.gene_name)

	else:
		s['gene_ids'] = '.'
		s['gene_names'] = '.'

	return s

cnv = cnv.apply(get_genes, axis = 1)
cnv.to_csv('call_tables/3_annotate_gencode.txt', sep = '\t', index = False)
