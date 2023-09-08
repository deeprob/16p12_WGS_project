
import pandas as pd

# I'll first annotate genes because the gene file is significantly smaller than the exon file
# After I've annotated and filtered for genes, I'll annotate exons

# Annotate gencode genes using gencode.v19.parsed.genes.csv

# gencode.v19.parsed.genes.csv was originally made by Anastasia
# Anastasia's parsing scripts are here: Dropbox/16p12.2 project/Human patients project/WGS paper/5_SNV pipelines-annotations/GENCODE annotations
genes = pd.read_csv('../../2022_02_14/gencode.v19.parsed.genes.csv')
genes = genes[genes.gene_type=='protein_coding']
print(genes)


# CNVs
cnv = pd.read_csv('bed_files/12_frequency_filter.bed', sep = '\t')

# Get genes
def get_genes(s):
	chr = s.Chr
	start = s.Start
	end = s.End

	# Get exons from gencode file that are at least partially contained within the CNV length
	gene_df = genes[(genes.Chrom==chr) & (((genes.Start>=start) & (genes.Start<=end)) | ((genes.End>=start) & (genes.End<=end)) | ((genes.Start<=start) & (genes.End>=end)))]

	if gene_df.shape[0]>0:
		# Get a list of gene IDs and gene names
		# NOTE: Some gene symbols have multiple IDs
		gene_ids = list(gene_df.gene_id.unique())
		s['gene_ids'] = ' '.join(gene_ids)
		gene_names = list(gene_df.gene_name.unique())
		s['gene_names'] = ' '.join(gene_names)

	else:
		s['gene_ids'] = '.'
		s['gene_names'] = '.'

	return s

cnv = cnv.apply(get_genes, axis = 1)
cnv = cnv[cnv.gene_ids != '.']
print(cnv)

# Annotate exons using gencode.v19.parsed.exons.csv
# gencode.v19.parsed.exons.csv was originally made by Anastasia
# Anastasia's parsing scripts are here: Dropbox/16p12.2 project/Human patients project/WGS paper/5_SNV pipelines-annotations/GENCODE annotations
exons = pd.read_csv('../../2022_02_14/gencode.v19.parsed.exons.csv')
exons = exons[exons.gene_type=='protein_coding']

# I decided to skip this step to clean up the output - we don't really need the exon numbers
# To speed up gathering attributes, limit exons to only those found in alread-annotated genes
#anno_genes = [lst.split(' ') for lst in cnv.gene_ids.to_list()]
#anno_lst = list(set([item for sublist in anno_genes for item in sublist]))
#exons = exons[exons.gene_id.isin(anno_lst)]

def get_attribute(row, attribute_wanted):
	attributes = row['Attribute']
	sattributes = attributes.split(';')
	for attribute in sattributes:
		attribute = attribute.strip()
		if attribute.startswith(attribute_wanted):
			if '"' in attribute:
				return attribute.split('"')[1]
			else:
				return attribute.split(' ')[1]
	return ''

#exons['exon_number'] = exons.apply(lambda row: get_attribute(row, 'exon_number'), axis = 1)
#exons['exon_id'] = exons.apply(lambda row: get_attribute(row, 'exon_id'), axis = 1)

def get_exons(s):
	chr = s.Chr
	start = s.Start
	end = s.End

	# Get exons from gencode file that are at least partially contained within the CNV length
	exon_df = exons[(exons.Chrom==chr) & (((exons.Start>=start) & (exons.Start<=end)) | ((exons.End>=start) & (exons.End<=end)) | ((exons.Start<=start) & (exons.End>=end)))]

	if exon_df.shape[0]>0:
		return True

	else:
		return False


cnv = cnv[cnv.apply(get_exons, axis = 1)]
print(cnv)

# Also, for downstream add a variant id
cnv['variant_id'] = cnv['Chr']+'_'+cnv['Start'].astype(str)+'_'+cnv['End'].astype(str)+'_'+cnv['Sample']

cnv.to_csv('bed_files/14_gene_filter.bed', sep = '\t', index = False)
