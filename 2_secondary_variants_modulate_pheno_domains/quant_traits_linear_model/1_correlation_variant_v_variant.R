


library(xlsx)         # to load excel sheets
library(effsize)	  # cohen's D
library(glue)


get_sample_type = function(s) {
	if (s == 'P') {
		return('Proband')
	}
	else if (s == 'MC') {
		return('Carrier Parent')
	}
	else if (s == 'FC') {
		return('Carrier Parent')
	}
	else if (s == 'FNC') {
		return('Noncarrier Parent')
	}
	else if (s == 'MNC') {
		return('Noncarrier Parent')
	}
	else {
		return('Other')
	}
}

dropbox_location = '~/Dropbox/'

filename = glue('{dropbox_location}/16p12.2 project/Human patients project/WGS paper/11_Variant Integration/16p12_cohort_summary_v17.xlsx')


df = read.xlsx(filename, sheetIndex=1, check.names=FALSE)
df['Relationship_Label'] = sapply(df[, 'Relationship'], get_sample_type)
rownames(df) = df$Sample

sample_groups = c('Proband', 'Carrier Parent', 'Noncarrier Parent')
variant_groups = c('Missense_CADD25', 'Missense_CADD25_LOEUF035', 'LOF', 'LOF_LOEUF_035', 'Splice_CADD25', 'Splice_CADD25_LOEUF035', 'genes_del', 'genes_dup', 'dels_loeuf', 'dups_loeuf', 'STRs_exonic', 'STRs_exonic_LOEUF035', 'Rare_Deleterious_SNVs', 'Rare_Deleterious_SNVs_LOEUF', 'enhancer', 'promoter', 'UTR5', 'intelligence_PRS', 'SCZ_PRS', 'educational_attainment_PRS', 'autism_PRS')

stats_df <- data.frame(Variant1=character(0), Variant2=character(0), Sample_Group=character(0), Pvalue=numeric(0), Num_samples=numeric(0), Coeff=numeric(0))
stats_df$Variant1 = as.character(stats_df$Variant1)
stats_df$Variant2 = as.character(stats_df$Variant2)
stats_df$Sample_Group = as.character(stats_df$Sample_Group)

for (sample_group in sample_groups) {
	subdf = df[df$Relationship_Label == sample_group,]


	for (variant_col1 in variant_groups){
		for (variant_col2 in variant_groups){

				tmp_mat = subdf[,c(variant_col1, variant_col2)]
				tmp_mat[,variant_col1] = as.numeric(tmp_mat[,variant_col1])
				tmp_mat[,variant_col2] = as.numeric(tmp_mat[,variant_col2])
				tmp_mat = na.omit(tmp_mat)
				
				num_samples = dim(tmp_mat)[1]

				
				a = tmp_mat[, variant_col1]
				b = tmp_mat[, variant_col2]
				res = cor.test(a, b)
				coeff = as.numeric(res$estimate)
				p_val = res$p.value
				

				app = c(variant_col1, variant_col2, sample_group, p_val, num_samples, coeff)
				stats_df[nrow(stats_df)+1,] = app


		}
	}
}



# FDR correction
stats_df['use_in_fdr_calculation'] = FALSE
for (variant_col1 in variant_groups) {
	for (variant_col2 in variant_groups) {
		if (variant_col1 == variant_col2) break
		
		stats_df[(stats_df$Variant1 == variant_col1) & (stats_df$Variant2 == variant_col2), 'use_in_fdr_calculation'] = T
	}
}



stats_df[,'FDR'] = NA
for (sample_group in sample_groups) {
	subdf = stats_df[stats_df$Sample_Group == sample_group,]
	subdf = subdf[subdf$use_in_fdr_calculation, ]
	
	pvalues = subdf$Pvalue
	adj_pvalues = p.adjust(pvalues, method='bonferroni')
	subdf[,'FDR'] = adj_pvalues
	
	for (index in rownames(subdf)){
		fdr = subdf[index, 'FDR']
		stats_df[index, 'FDR'] = fdr
		
		variant_col1 = subdf[index, 'Variant1']
		variant_col2 = subdf[index, 'Variant2']
		
		stats_df[(stats_df$Variant2 == variant_col1) &
			(stats_df$Variant1 == variant_col2) &
			(stats_df$Sample_Group == sample_group)
			, 'FDR'] = fdr
		
		
		
	}
}
		





write.csv(stats_df, 'statistics/variant_v_variant.csv', row.names=F)

























