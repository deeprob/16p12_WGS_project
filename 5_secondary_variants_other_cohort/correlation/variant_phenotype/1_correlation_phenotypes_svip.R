



library(effsize)	  # cohen's D
library(glue)


cohorts = c('16p-deletion', '16p-duplication')
variant_cols = c('Missense CADD>25', 'CADD25 LOEUF<0.35', 'LOF CADD25_NA', 'LOF CADD LOEUF<0.35', 'Splice CADD25', 'Splice CADD25 LOEUF<0.35', 'Rare_Deleterious_SNVs', 'Rare_Deleterious_SNVs_LOEUF', 'Genes DEL', 'Genes DUP', 'DELs LOEUF<0.35', 'DUPs LOEUF<0.35', 'SCZ_PRS', 'intelligence_PRS', 'educational_attainment_PRS', 'autism_PRS')
pheno_cols = c('Full-scale IQ','ABCL_CBCL external','ABCL_CBCL internal','SRS raw','BSI','BMI Z-score','Head circumference Z-score')


stats_df <- data.frame(Variant=character(0), Phenotype=character(0), First_hit=character(0), Pvalue=numeric(0), Num_samples=numeric(0), Coeff=numeric(0))
stats_df$Variant = as.character(stats_df$Variant)
stats_df$Phenotype = as.character(stats_df$Phenotype)
stats_df$First_hit = as.character(stats_df$First_hit)

for (cohort in cohorts) {
	filename = glue('../1_variant_preperation/SVIP.csv')
	df = read.csv(filename, check.names=F)
	df = df[df[,'Family type'] == cohort,]
	colnames(df) = gsub('/', '_', colnames(df))
	rownames(df) = df$Sample
	


	for (variant_col in variant_cols){
		for (pheno_col in pheno_cols) {

			tmp_mat = df[,c(variant_col, pheno_col)]
			tmp_mat[,variant_col] = as.numeric(tmp_mat[,variant_col])
			tmp_mat[,pheno_col] = as.numeric(tmp_mat[,pheno_col])
			tmp_mat = na.omit(tmp_mat)
			
			num_samples = dim(tmp_mat)[1]

			
			a = tmp_mat[, variant_col]
			b = tmp_mat[, pheno_col]
			res = cor.test(a, b)
			coeff = res$estimate
			p_val = res$p.value
			

			app = c(variant_col, pheno_col, cohort, p_val, num_samples, coeff)
			stats_df[nrow(stats_df)+1,] = app
		}

	}
}



# FDR correction
stats_df[,'FDR'] = 1.0
for (pheno_col in pheno_cols) {
	for (cohort in cohorts) {
		subdf = stats_df[stats_df$First_hit == cohort,]
		subdf = subdf[subdf$Phenotype == pheno_col,]
		
		pvalues = subdf$Pvalue
		adj_pvalues = p.adjust(pvalues, method='bonferroni')
		subdf[,'FDR'] = adj_pvalues
		
		for (index in rownames(subdf)){
			fdr = subdf[index, 'FDR']
			stats_df[index, 'FDR'] = fdr
		}
	}
}

colnames(stats_df) = c("Variant class/PRS", "Phenotype", "Group", "P-value", "Sample size", "Direction", "FDR")

print(head(stats_df))

write.csv(stats_df, 'statistics/variant_v_phenotype.csv', row.names=F)



















