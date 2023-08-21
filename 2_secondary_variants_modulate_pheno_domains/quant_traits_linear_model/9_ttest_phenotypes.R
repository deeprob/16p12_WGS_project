


library(xlsx)         # to load excel sheets
library(effsize)	  # cohen's D
library(glue)
library(sjstats)      # partial omega squared test




dropbox_location = '~/Dropbox/'

filename = glue('{dropbox_location}/16p12.2 project/Human patients project/WGS paper/11_Variant Integration/16p12_cohort_summary_v17.xlsx')


df = read.xlsx(filename, sheetIndex=1, check.names=FALSE)
rownames(df) = df$Sample

# FILTER
df = df[df$Relationship == 'P',]

# for each phenotype split the cohort in two
# set to 0, 1, or N/A
numerical2binary = function(a, split_value=0) {
	if (is.na(a)) {
		return(NA)
	} else if (a<=split_value) {
		return(0)
	} else if (a>split_value) {
		return(1)
	} 
}

col = 'Child_ID_DD'
df[,col] = unlist(sapply(df[, col], numerical2binary, split_value=2))

col = 'Child_behav'
df[,col] = unlist(sapply(df[, col], numerical2binary, split_value=1))

col = 'Child_psych'
df[,col] = unlist(sapply(df[, col], numerical2binary, split_value=0))

col = 'Child_nervous_system'
df[,col] = unlist(sapply(df[, col], numerical2binary, split_value=0))

col = 'Child_congenital'
df[,col] = unlist(sapply(df[, col], numerical2binary, split_value=1))

col = 'Child_craniofacial'
df[,col] = unlist(sapply(df[, col], numerical2binary, split_value=2))

col = 'De_vrie'
df[,col] = unlist(sapply(df[, col], numerical2binary, split_value=7))




phenotype_groups = c('Child_ID_DD', 'Child_behav', 'Child_psych', 'Child_nervous_system', 'Child_congenital', 'Child_craniofacial', 'De_vrie')
variant_groups = c('Missense_CADD25', 'Missense_CADD25_LOEUF035', 'LOF', 'LOF_LOEUF_035', 'Splice_CADD25', 'Splice_CADD25_LOEUF035', 'genes_del', 'genes_dup', 'dels_loeuf', 'dups_loeuf', 'STRs_exonic', 'STRs_exonic_LOEUF035', 'Rare_Deleterious_SNVs', 'Rare_Deleterious_SNVs_LOEUF', 'enhancer', 'promoter', 'UTR5', 'intelligence_PRS', 'SCZ_PRS', 'educational_attainment_PRS', 'autism_PRS')


stats_df <- data.frame(Variant=character(0), Phenotype=character(0), Pvalue=numeric(0), Num_samples=numeric(0), Cohen_D=numeric(0))
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

		
		formula = as.formula(glue('`{variant_col}` ~ `{phenotype_col}`'))
		res = t.test(formula, data=tmp_mat)
		pval = res$p.value
		
		a = tmp_mat[tmp_mat[,phenotype_col] == 1,variant_col]
		b = tmp_mat[tmp_mat[,phenotype_col] == 0,variant_col]
						
		res = cohen.d(a, b)			
		cohen_d    = res$estimate
		

		app = c(variant_col, phenotype_col, pval, num_samples, cohen_d)
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





write.csv(stats_df, 'statistics/ttest_phenotypes.csv', row.names=F)

























