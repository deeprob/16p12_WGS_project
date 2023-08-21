#!/bin/R


library(xlsx)         # to load excel sheets
library(glue)         # to format strings



dropbox_location = '~/Dropbox'




# load in master table
filename = glue('{dropbox_location}/16p12.2 project/Human patients project/WGS paper/11_Variant Integration/16p12_cohort_summary_v17.xlsx')

df = read.xlsx(filename, sheetIndex=1, check.names=FALSE)

# rename rows
rownames(df) = df$Sample

# FILTER 1
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





pheno_cols = c('Child_ID_DD', 'Child_behav', 'Child_psych', 'Child_nervous_system', 'Child_congenital', 'Child_craniofacial', 'De_vrie')
prs_cols = c('SCZ_PRS', 'intelligence_PRS', 'educational_attainment_PRS', 'autism_PRS')



df['All_rare_del_var'] = df[,'Missense_CADD25'] + df[,'LOF'] +
	df[,'Splice_CADD25'] + df[, 'genes_del'] + df[, 'genes_dup'] +
	df[,'STRs_exonic']

# create a dataframe to populate with statistics
summ_df = data.frame(matrix(ncol=6, nrow=0))
colnames(summ_df) = c('Phenotype', 'PRS', 'pvalue', 'FDR', 'R', 'num_samples')
summ_df = summ_df[-c(1),]



for (pheno_col in pheno_cols) {
	
		for (prs_col in prs_cols) {


			tmp_mat = df[,c(pheno_col, prs_col, 'All_rare_del_var')]
			tmp_mat = na.omit(tmp_mat)
			tmp_mat = tmp_mat[tmp_mat[,pheno_col] == 1,]
			
			num_samples = dim(tmp_mat)[1]
			
			a = tmp_mat[, prs_col]
			b = tmp_mat[, 'All_rare_del_var']
			
			res = cor.test(a, b)
			p_val = res$p.value
			coeff = as.numeric(res$estimate)
			
			
			append = c(pheno_col, prs_col, p_val, 0, coeff, num_samples)
			print(append)
			
			
			# summ_df = rbind(summ_df, append)
			summ_df[nrow(summ_df) + 1,] = append
		}

}




# FDR correction
pvalues = summ_df[,'pvalue']
summ_df[,'FDR'] = p.adjust(pvalues, method='bonferroni')


# save statistics
filename = 'statistics/correlations_rare_variant_prs.tsv'
write.table(summ_df, filename, sep='\t', row.names=F)



