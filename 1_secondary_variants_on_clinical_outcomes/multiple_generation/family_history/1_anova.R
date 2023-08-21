


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


phenotype_groups = c('max_clusters_new','sum_clusters_new','mother_clusters_new','father_clusters_new','Carrier_clusters_new','Noncarrier_clusters_new')
variant_groups = c('Missense_CADD25', 'Missense_CADD25_LOEUF035', 'LOF', 'LOF_LOEUF_035', 'Splice_CADD25', 'Splice_CADD25_LOEUF035', 'genes_del', 'genes_dup', 'dels_loeuf', 'dups_loeuf', 'STRs_exonic', 'STRs_exonic_LOEUF035', 'Rare_Deleterious_SNVs', 'Rare_Deleterious_SNVs_LOEUF', 'enhancer', 'promoter', 'UTR5', 'intelligence_PRS', 'SCZ_PRS', 'educational_attainment_PRS', 'autism_PRS')


stats_df <- data.frame(Variant=character(0), Family_History=character(0), Pvalue=numeric(0), Num_samples=numeric(0), Partial_Omega_Sq=numeric(0))
stats_df$Variant = as.character(stats_df$Variant)
stats_df$Family_History = as.character(stats_df$Phenotype)

for (phenotype_col in phenotype_groups) {
	subdf = df[!is.na(df[,phenotype_col]),]
	subdf = subdf[!(subdf[,phenotype_col]==''),]

	for (variant_col in variant_groups){
		tmp_mat = subdf[,c(variant_col, phenotype_col)]
		tmp_mat[,variant_col] = as.numeric(tmp_mat[,variant_col])
		tmp_mat[,phenotype_col] = as.factor(tmp_mat[,phenotype_col])
		tmp_mat = na.omit(tmp_mat)
		
		num_samples = dim(tmp_mat)[1]

		
		formula = as.formula(glue('`{variant_col}` ~ `{phenotype_col}`'))
		res = aov(as.formula(formula), tmp_mat)
		pval = summary(res)[[1]][["Pr(>F)"]][1]
		
		o_res = anova_stats(res)
		partial_om = o_res[,'partial.omegasq'][1]
		

		app = c(variant_col, phenotype_col, pval, num_samples, partial_om)
		stats_df[nrow(stats_df)+1,] = app


	}
}


print(head(stats_df))
print(tail(stats_df))


# FDR correction
stats_df[,'FDR'] = NA
for (phenotype_group in phenotype_groups) {
	
	# do FDR correction once for each phenotype
	subdf = stats_df[stats_df$Family_History == phenotype_group,]
	
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





write.csv(stats_df, 'statistics/anova_family_history.csv', row.names=F)

























