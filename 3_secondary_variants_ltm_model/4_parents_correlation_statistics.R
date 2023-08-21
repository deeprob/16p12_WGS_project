#!/bin/R


library(xlsx)         # to load excel sheets
library(glue)         # to format strings



dropbox_location = '~/Dropbox'






get_sample_type = function(s) {
	if (s == 'P') {
		return('Proband')
	}
	else if (s == 'MC') {
		return('Carrier_parent')
	}
	else if (s == 'FC') {
		return('Carrier_parent')
	}
	else if (s == 'FNC') {
		return('Noncarrier_parent')
	}
	else if (s == 'MNC') {
		return('Noncarrier_parent')
	}
	else {
		return('Other')
	}
}



# load in master table
filename = glue('{dropbox_location}/16p12.2 project/Human patients project/WGS paper/11_Variant Integration/16p12_cohort_summary_v17.xlsx')
df = read.xlsx(filename, sheetIndex=1, check.names=FALSE)

# rename rows
rownames(df) = df$Sample


df['Relationship_Label'] = sapply(df[, 'Relationship'], get_sample_type)
df = df[df$Relationship_Label %in% c('Proband', 'Carrier_parent', 'Noncarrier_parent'),]


prs_cols = c('SCZ_PRS', 'intelligence_PRS', 'educational_attainment_PRS', 'autism_PRS')



df['All_rare_del_var'] = df[,'Missense_CADD25'] + df[,'LOF'] +
	df[,'Splice_CADD25'] + df[, 'genes_del'] + df[, 'genes_dup'] +
	df[,'STRs_exonic']

# create a dataframe to populate with statistics
summ_df = data.frame(matrix(ncol=6, nrow=0))
colnames(summ_df) = c('Role', 'PRS', 'pvalue', 'FDR', 'R', 'num_samples')
summ_df = summ_df[-c(1),]


for (relationship in c('Proband', 'Carrier_parent', 'Noncarrier_parent')) {
	print(relationship)
	
	for (prs_col in prs_cols) {
		tmp_mat = df[,c('Relationship_Label', prs_col, 'All_rare_del_var')]
		tmp_mat = na.omit(tmp_mat)
		tmp_mat = tmp_mat[tmp_mat[,'Relationship_Label'] == relationship,]
		print(head(tmp_mat))
		
		num_samples = dim(tmp_mat)[1]
		
		a = tmp_mat[, prs_col]
		b = tmp_mat[, 'All_rare_del_var']
		
		res = cor.test(a, b)
		p_val = res$p.value
		coeff = as.numeric(res$estimate)
		
		
		append = c(relationship, prs_col, p_val, 0, coeff, num_samples)
		
		
		
		summ_df[nrow(summ_df) + 1,] = append

	}

}


summ_df[,'FDR'] = NA
for (relationship in c('Proband', 'Carrier_parent', 'Noncarrier_parent')) {
	subdf = summ_df[summ_df$Role == relationship,]
	
	pvalues = subdf$pvalue
	adj_pvalues = p.adjust(pvalues, method='bonferroni')
	subdf[,'FDR'] = adj_pvalues
	
	for (index in rownames(subdf)){
		fdr = subdf[index, 'FDR']
		summ_df[index, 'FDR'] = fdr
	}
}


# FDR correction
# pvalues = summ_df[,'pvalue']
# summ_df[,'FDR'] = p.adjust(pvalues, method='bonferroni')


# save statistics
filename = 'statistics/compared_to_parents_correlations.tsv'
write.table(summ_df, filename, sep='\t', row.names=F)



