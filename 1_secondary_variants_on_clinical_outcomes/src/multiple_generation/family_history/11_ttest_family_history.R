


library(xlsx)         # to load excel sheets
library(effsize)	  # cohen's D
library(glue)
library(sjstats)      # partial omega squared test




dropbox_location = 'D:/Dropbox/16p12.2 project'

filename = glue('{dropbox_location}/Human patients project/WGS paper/11_Variant Integration/16p12_cohort_summary_v21.csv')


df = read.csv(filename, check.names=FALSE)
rownames(df) = df$Sample

# FILTER
df = df[df$Relationship == 'P',]

# Add in additional family history clusters
fam_hist2 <- read.csv(glue('{dropbox_location}/Human patients project/WGS paper/37_Family History/Analysis_files/3_grouped_clusters.csv'), check.names = F)
fam_hist2 <- fam_hist2[, c('Proband', "number_parents", "phenotype_parents")]
colnames(fam_hist2) <- c('Sample', "number_parents", "phenotype_parents")
df <- merge(df, fam_hist2, by = 'Sample', all.x=T, all.y=F)


#----------------------------
# Prepare phenotypes
#----------------------------

family_history_cols = c('NEWmax_clusters','NEWsum_clusters','NEWmother_clusters','NEWfather_clusters','NEWCarrier_clusters','NEWNoncarrier_clusters', 'NEWmother_father_clusters', 'NEWCarrier_Noncarrier_clusters', "number_parents", "phenotype_parents")
family_history_groups = list(c("A","B"), c("A","C"), c("B","C"))

new_family_history_groups = function (s, group1, group2){
	if (is.na(s)) {
		return(NA)
	} else if (s == group1) {
		return(0)
	} else if (s == group2) {
		return(1)
	} 
	return(NA)
}

for (col in family_history_cols) {

	
	for (group in family_history_groups) {
		suffix = glue_collapse(group, sep='_v_')
		new_col = glue('{col}_{suffix}')
		
		group1 = group[1]
		group2 = group[2]
		
		df[,new_col] = unlist(sapply(df[, col], new_family_history_groups, group1=group1, group2=group2))


		
	}
}


phenotype_groups <- colnames(df)[64:93]
# phenotype_groups = c('max_clusters_new_A_v_B','max_clusters_new_A_v_C','max_clusters_new_B_v_C','sum_clusters_new_A_v_B','sum_clusters_new_A_v_C','sum_clusters_new_B_v_C','mother_clusters_new_A_v_B','mother_clusters_new_A_v_C','mother_clusters_new_B_v_C','father_clusters_new_A_v_B','father_clusters_new_A_v_C','father_clusters_new_B_v_C','Carrier_clusters_new_A_v_B','Carrier_clusters_new_A_v_C','Carrier_clusters_new_B_v_C','Noncarrier_clusters_new_A_v_B','Noncarrier_clusters_new_A_v_C','Noncarrier_clusters_new_B_v_C')
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





write.csv(stats_df, 'statistics/ttest_family_history.csv', row.names=F)

























