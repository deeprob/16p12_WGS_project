



library(effsize)	  # cohen's D
library(glue)


cohorts = c('16p-deletion', '16p-duplication')
variant_cols = c('Missense CADD>25', 'CADD25 LOEUF<0.35', 'LOF CADD25_NA', 'LOF CADD LOEUF<0.35', 'Splice CADD25', 'Splice CADD25 LOEUF<0.35', 'Rare_Deleterious_SNVs', 'Rare_Deleterious_SNVs_LOEUF', 'Genes DEL', 'Genes DUP', 'DELs LOEUF<0.35', 'DUPs LOEUF<0.35', 'SCZ_PRS', 'intelligence_PRS', 'educational_attainment_PRS', 'autism_PRS')




stats_df <- data.frame(Variant1=character(0), Variant2=character(0), First_hit=character(0), Pvalue=numeric(0), Num_samples=numeric(0), Coeff=numeric(0))
stats_df$Variant1 = as.character(stats_df$Variant1)
stats_df$Variant2 = as.character(stats_df$Variant2)
stats_df$First_hit = as.character(stats_df$First_hit)


for (cohort in cohorts) {
	filename = glue('../1_variant_preperation/SVIP.csv')
	df = read.csv(filename, check.names=F)
	df = df[df[,'Family type'] == cohort,]
	colnames(df) = gsub('/', '_', colnames(df))
	rownames(df) = df$Sample


	for (variant_col1 in variant_cols){
		for (variant_col2 in variant_cols){

				tmp_mat = df[,c(variant_col1, variant_col2)]
				tmp_mat[,variant_col1] = as.numeric(tmp_mat[,variant_col1])
				tmp_mat[,variant_col2] = as.numeric(tmp_mat[,variant_col2])
				tmp_mat = na.omit(tmp_mat)
				
				num_samples = dim(tmp_mat)[1]

				
				a = tmp_mat[, variant_col1]
				b = tmp_mat[, variant_col2]
				res = cor.test(a, b)
				coeff = res$estimate
				p_val = res$p.value
				

				app = c(variant_col1, variant_col2, cohort, p_val, num_samples, coeff)
				stats_df[nrow(stats_df)+1,] = app


		}
	}
}


# FDR correction
stats_df['use_in_fdr_calculation'] = FALSE
for (variant_col1 in variant_cols) {
	for (variant_col2 in variant_cols) {
		if (variant_col1 == variant_col2) break
		
		stats_df[(stats_df$Variant1 == variant_col1) & (stats_df$Variant2 == variant_col2), 'use_in_fdr_calculation'] = T
	}
}


stats_df[,'FDR'] = NA
for (group in cohorts) {
	subdf = stats_df[stats_df$First_hit == group,]
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
			(stats_df$First_hit == group)
			, 'FDR'] = fdr
		
		
		
	}
}
		





write.csv(stats_df, 'statistics/variant_v_variant.csv', row.names=F)

























