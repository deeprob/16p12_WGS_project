


library(xlsx)         # to load excel sheets
library(effsize)	  # cohen's D
library(glue)
library(sjstats)      # partial omega squared test




dropbox_location = 'C:/Users/mattq/Dropbox'

filename = glue('{dropbox_location}/Human patients project/WGS paper/11_Variant Integration/16p12_cohort_summary_v21.xlsx')


df = read.xlsx(filename, sheetIndex=1, check.names=FALSE)
rownames(df) = df$Sample

# FILTER
df = df[df$Relationship == 'P',]


phenotype_groups = c('bmi_z_score','head_circumference_z_score','hrs_mat','srs')
variant_groups = c('Missense_CADD25', 'Missense_CADD25_LOEUF035', 'LOF', 'LOF_LOEUF_035', 'Splice_CADD25', 'Splice_CADD25_LOEUF035', 'genes_del', 'genes_dup', 'dels_loeuf', 'dups_loeuf', 'STRs_exonic', 'STRs_exonic_LOEUF035', 'Rare_Deleterious_SNVs', 'Rare_Deleterious_SNVs_LOEUF', 'enhancer', 'promoter', 'UTR5', 'intelligence_PRS', 'SCZ_PRS', 'educational_attainment_PRS', 'autism_PRS')


stats_df <- data.frame(Variant=character(0), Phenotype=character(0), Pvalue=numeric(0), Num_samples=numeric(0), Coeff=numeric(0))
stats_df$Variant = as.character(stats_df$Variant)
stats_df$Phenotype = as.character(stats_df$Phenotype)


for (phenotype_col in phenotype_groups) {
	subdf = df[!is.na(df[,phenotype_col]),]

	for (variant_col in variant_groups){
		tmp_mat = subdf[,c(variant_col, phenotype_col)]
		tmp_mat[,variant_col] = as.numeric(tmp_mat[,variant_col])
		tmp_mat[,phenotype_col] = as.numeric(tmp_mat[,phenotype_col])
		tmp_mat = na.omit(tmp_mat)
		
		num_samples = dim(tmp_mat)[1]

		
		a = tmp_mat[, variant_col]
		b = tmp_mat[, phenotype_col]
		res = cor.test(a, b)
		coeff = as.numeric(res$estimate)
		p_val = res$p.value
		

		app = c(variant_col, phenotype_col, p_val, num_samples, coeff)
		stats_df[nrow(stats_df)+1,] = app


	}
}


print(head(stats_df))
print(tail(stats_df))


# FDR correction
stats_df[,'FDR'] = NA
for (phenotype_group in phenotype_groups) {
	
	# do FDR correction once for each phenotype
	subdf = stats_df[stats_df$Phenotype == phenotype_group,]
	
	pvalues = subdf$Pvalue
	adj_pvalues = p.adjust(pvalues, method='bonferroni')
	subdf[,'FDR'] = adj_pvalues
	
	for (index in rownames(subdf)){
		fdr = subdf[index, 'FDR']
		stats_df[index, 'FDR'] = fdr
		
	}
}
		


print(head(stats_df))
print(tail(stats_df))





write.csv(stats_df, 'statistics/correlation_phenotypes.csv', row.names=F)

























