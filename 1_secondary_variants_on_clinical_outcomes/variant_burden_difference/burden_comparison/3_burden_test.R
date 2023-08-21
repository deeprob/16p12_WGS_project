
library(xlsx)         # to load excel sheets
library(glue)         # to format strings
library(stringr)      # for str_replace_all
library(effsize)	  # cohen's D

#----------------------------
# helper functions
#----------------------------
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




#----------------------------
# load in data
#----------------------------



dropbox_location = 'C:/Users/corny/Dropbox/'

filename = glue('{dropbox_location}/16p12.2 project/Human patients project/WGS paper/11_Variant Integration/16p12_cohort_summary_v17.xlsx')

df = read.xlsx(filename, sheetIndex=1, check.names=FALSE)


rownames(df) = df$Sample



# FILTER 1
df['Relationship_Label'] = sapply(df[, 'Relationship'], get_sample_type)
df = df[df$Relationship_Label %in% c('Proband', 'Carrier Parent'),]


filename = glue('{dropbox_location}/16p12.2 project/Human patients project/WGS paper/11_Variant Integration/16p12_cohort_summary_v17.xlsx')



estonia_df = read.xlsx(glue('{dropbox_location}/16p12.2 project/Human patients project/WGS paper/11_Variant Integration/Estonian_Summary.xlsx'), sheetIndex=1, check.names=FALSE)
# FILTER
estonia_df = estonia_df[estonia_df$Relationship %in% c('P'),]



#---------------------------------
# do burden tests
#---------------------------------


# groups to test
relationships = c('Proband', 'Carrier Parent')
variant_groups = c('Missense_CADD25', 'Missense_CADD25_LOEUF035', 'LOF', 'LOF_LOEUF_035', 'Splice_CADD25', 'Splice_CADD25_LOEUF035', 'genes_del', 'genes_dup', 'dels_loeuf', 'dups_loeuf', 'STRs_exonic', 'STRs_exonic_LOEUF035', 'Rare_Deleterious_SNVs', 'Rare_Deleterious_SNVs_LOEUF', 'enhancer', 'promoter', 'UTR5', 'intelligence_PRS', 'SCZ_PRS', 'educational_attainment_PRS', 'autism_PRS')




# store results here
cols = c('P-value','statistic','group1','Comparison','group2','Variant Group', 'Cohens D', 'Sample_Estonia', 'Sample_group')
summ = data.frame(matrix(ncol=length(cols)))
colnames(summ) = cols
summ = summ[-c(1),]



# for each variant group calculate t-test between each group
for (variant_group in variant_groups) {
	for (group in relationships) {
			
			
		# FILTER
		subdf = df[!is.na(df[,variant_group]),]
		sub_estonia_df = estonia_df[!is.na(estonia_df[,variant_group]),]
		
		# FILTER
		subdf = subdf[subdf$Relationship_Label == group,]
		
		a = sub_estonia_df[,variant_group]
		b = subdf[,variant_group]
		
		res = wilcox.test(a, b)
		
		pval  = res$p.value
		tstat = as.numeric(res$statistic)
		
		res2 = cohen.d(a, b, na.rm=TRUE)			
		cohen_d  = res2$estimate
		
		summ[nrow(summ) + 1,] = c(pval, tstat, 'Estonia', 'different than', group, variant_group, cohen_d, length(a), length(b))

	}
}

#---------------------------------
# Filter statistics and bonferroni correction
#---------------------------------


# FDR correction
summ[,'FDR'] = p.adjust(summ$`P-value`, method='bonferroni')


# file to write out statistics
outfile = 'statistics/wilcox_test.csv'

# save
write.csv(summ, file=outfile, row.names=F)

