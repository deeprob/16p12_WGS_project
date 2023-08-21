#!/bin/python3





library(glue)
library(dplyr)


dropbox = '~/Dropbox'


# load in scrnaseq data
scrna_df = read.csv('allen_institute_tables/cell_type_specific_expression_subclass.csv')
scrna_df = rename(scrna_df, feature = X)
rownames(scrna_df) = scrna_df$feature
cell_types = colnames(scrna_df)
cell_types = cell_types[2:length(cell_types)]

# load in genes with second hits
variant_classes = c('rare_deleterious_snvs', 'dups', 'dels', 'strs_std_2')
columns = c('Sample', 'Gene', 'variant_class')
variant_df = data.frame(matrix(ncol = length(columns), nrow = 0))
colnames(variant_df) = columns

for (variant_class in variant_classes) {
	app = read.csv(glue('../../19_Variant_GO_enrichment/figure_for_16p12_cohort/variants/{variant_class}.csv'))
	app[,'variant_class'] = variant_class
	variant_df = rbind(variant_df, app)
}


# restrict to probands
master_df = read.csv('../../11_Variant Integration/16p12_cohort_summary_v21.csv')
# master_df = master_df.set_index('Sample', drop=False)
# only want samples with WGS 
master_df = master_df[!is.na(master_df$Missense_CADD25),]
master_df = master_df[master_df$Relationship == 'P',]
variant_df = variant_df[variant_df$Sample %in% master_df$Sample,]


# use the intersection of genes in GENCODE and 
# genes in the scrnaseq table as the gene universe
gencode_df = read.csv('../../5_SNV pipelines-annotations/GENCODE annotations/gencode.v19.parsed.genes.csv')
gencode_df = gencode_df[gencode_df$gene_type == 'protein_coding',]
gencode_genes = unique(gencode_df$gene_name)
scrna_genes = unique(rownames(scrna_df))
gene_universe = intersect(gencode_genes, scrna_genes)


# restrict scrna_df and variant_df to gene universe
scrna_df = scrna_df[gene_universe,]
variant_df = variant_df[variant_df$Gene %in% gene_universe,]


# for each variant class do a fishers exact test for each cell type
columns=c('variant_class', 'cell_type', 'oddsratio', 'pvalue', 'conf_int_lower', 'conf_int_upper', 'second_hit_and_pref', 'second_hit_and_not_pref', 'not_second_hit_and_pref', 'not_second_hit_and_pref')
stats = data.frame(matrix(ncol=length(columns), nrow=0))
colnames(stats) = columns
for (variant_class in variant_classes) {
	sub_variant_df = variant_df[variant_df$variant_class == variant_class,]
	for (cell_type in cell_types) {
		genes_with_second_hits = sub_variant_df$Gene
		genes_for_cell_type = rownames(scrna_df[scrna_df[,cell_type] == 1,])
		
		second_hit_and_pref = unique(intersect(genes_with_second_hits, genes_for_cell_type))
		second_hit_and_not_pref = unique(genes_with_second_hits[!(genes_with_second_hits %in% genes_for_cell_type)])
		not_second_hit_and_pref = unique(genes_for_cell_type[!(genes_for_cell_type %in% genes_with_second_hits)])
		not_second_hit_and_not_pref_tmp = unique(gene_universe[!(gene_universe %in% genes_with_second_hits)])
		not_second_hit_and_not_pref = unique(not_second_hit_and_not_pref_tmp[!(not_second_hit_and_not_pref_tmp %in% genes_for_cell_type)])
		print(length(second_hit_and_pref)+length(second_hit_and_not_pref)+length(not_second_hit_and_pref)+length(not_second_hit_and_not_pref))
		
		contingency_table = matrix(nrow=2, ncol=2)
		contingency_table[1,1] = length(second_hit_and_pref)
		contingency_table[1,2] = length(second_hit_and_not_pref)
		contingency_table[2,1] = length(not_second_hit_and_pref)
		contingency_table[2,2] = length(not_second_hit_and_not_pref)
		
		result = fisher.test(contingency_table)
		oddsratio = result$estimate
		pvalue = result$p.value
		conf_int_lower = result$conf.int[1]
		conf_int_upper = result$conf.int[2]
		
		stats[nrow(stats) + 1,] = c(variant_class, cell_type, oddsratio, pvalue, conf_int_lower, conf_int_upper, length(second_hit_and_pref), length(second_hit_and_not_pref), length(not_second_hit_and_pref), length(not_second_hit_and_not_pref))
		
	}
}



# FDR correction for each variant class (BH correction)
stats[,'FDR'] = NA
for (variant_class in variant_classes) {
	sub_stats = stats[stats$variant_class == variant_class,]
	sub_stats$FDR = p.adjust(sub_stats$pvalue, method='BH')
	for (i in rownames(sub_stats)) {
		stats[i, 'FDR'] = sub_stats[i, 'FDR']
	}
}




write.csv(stats, 'statistics/fishers_exact_with_R_subclass.csv', row.names=F)

















