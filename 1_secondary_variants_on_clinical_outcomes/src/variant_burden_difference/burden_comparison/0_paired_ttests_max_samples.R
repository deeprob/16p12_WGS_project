
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



dropbox_location = '~/Dropbox/'

filename = '../../11_Variant Integration/16p12_cohort_summary_v21.xlsx'

df = read.xlsx(filename, sheetIndex=1, check.names=FALSE)


rownames(df) = df$Sample



# FILTER 1
df['Relationship_Label'] = sapply(df[, 'Relationship'], get_sample_type)
df = df[df$Relationship_Label != 'Other',]

# FILTER
# remove samples if two non carrier parents
families_of_non_carrier_parents = df[df$Relationship_Label == 'Noncarrier Parent',]$Family
families_with_two_noncarrier_parents = families_of_non_carrier_parents[duplicated(families_of_non_carrier_parents)]
df = df[!(df$Family %in% families_with_two_noncarrier_parents & df$Relationship_Label == 'Noncarrier Parent'),]

# sort by family
df = df[order(df$Family), ]


#---------------------------------
# do t-tests
#---------------------------------


# groups to test
relationships = c('Proband', 'Carrier Parent', 'Noncarrier Parent')
variant_groups = c('UTR5', 'promoter', 'enhancer', 'Rare_Deleterious_SNVs_LOEUF', 'Rare_Deleterious_SNVs', 'STRs_exonic_LOEUF035', 'STRs_exonic', 'dups_loeuf', 'dels_loeuf', 'genes_dup', 'genes_del', 'Splice_CADD25_LOEUF035', 'Splice_CADD25', 'LOF_LOEUF_035', 'LOF', 'Missense_CADD25_LOEUF035', 'Missense_CADD25')
prs_groups = c('autism_PRS', 'educational_attainment_PRS', 'SCZ_PRS', 'intelligence_PRS')


# store results here
cols_rv = c('P-value','t-statistic','Conf. Int. 95% lower','Conf. Int. 95% upper','group1','Comparison','group2','Variant Group', 'Cohens D', 'Sample_size')
summ_rv = data.frame(matrix(ncol=length(cols_rv)))
colnames(summ_rv) = cols_rv
summ_rv = summ_rv[-c(1),]



# for each variant group calculate t-test between each group
for (variant_group in variant_groups) {
	for (group1 in relationships) {
		for(group2 in relationships) {
			
			if (group1 == group2) next
			
			# FILTER
			subdf = df[!is.na(df[,variant_group]),]
			
			# FILTER
			families_group1 = subdf[subdf$Relationship_Label == group1,]$Family
			families_group2 = subdf[subdf$Relationship_Label == group2,]$Family
			families = intersect(families_group1, families_group2)
			subdf = subdf[subdf$Family %in% families,]
			
			a = subdf[subdf$Relationship_Label == group1,variant_group]
			b = subdf[subdf$Relationship_Label == group2,variant_group]
			
			num_samples = length(a) + length(b)
			
			res = t.test(a,b, paired=TRUE, alternative='greater')
			
			pval  = res$p.value
			tstat = res$statistic
			conf  = res$conf.int[1]
			
			res2 = cohen.d(a, b, paired=TRUE, na.rm=TRUE)			
			cohen_d    = res2$estimate
			
			summ_rv[nrow(summ_rv) + 1,] = c(pval, tstat, conf, '', group1, 'greater than', group2, variant_group, cohen_d, num_samples)
			
		}
	}	
}

#---------------------------------
# PRS uses a two sided test and has different filters
#---------------------------------

dropbox_location = '~/Dropbox/'

filename = '../..//11_Variant Integration/16p12_cohort_summary_v21.xlsx'

df = read.xlsx(filename, sheetIndex=1, check.names=FALSE)


rownames(df) = df$Sample



# FILTER 1
df['Relationship_Label'] = sapply(df[, 'Relationship'], get_sample_type)
df = df[df$Relationship_Label != 'Other',]

# FILTER
# remove samples if two non carrier parents
families_of_non_carrier_parents = df[df$Relationship_Label == 'Noncarrier Parent',]$Family
families_with_two_noncarrier_parents = families_of_non_carrier_parents[duplicated(families_of_non_carrier_parents)]
df = df[!(df$Family %in% families_with_two_noncarrier_parents & df$Relationship_Label == 'Noncarrier Parent'),]

# sort by family
df = df[order(df$Family), ]

# store results here
cols_prs = c('P-value','t-statistic','Conf. Int. 95% lower','Conf. Int. 95% upper','group1','Comparison','group2','Variant Group', 'Cohens D', 'Sample_size')
summ_prs = data.frame(matrix(ncol=length(cols_prs)))
colnames(summ_prs) = cols_prs
summ_prs = summ_prs[-c(1),]


# for each variant group calculate t-test between each group
for (variant_group in prs_groups) {
	for (group1 in relationships) {
		for(group2 in relationships) {
			
			if (group1 == group2) next
			
			# FILTER
			subdf = df[!is.na(df[,variant_group]),]
			
			# FILTER
			families_group1 = subdf[subdf$Relationship_Label == group1,]$Family
			families_group2 = subdf[subdf$Relationship_Label == group2,]$Family
			families = intersect(families_group1, families_group2)
			subdf = subdf[subdf$Family %in% families,]
			
			a = subdf[subdf$Relationship_Label == group1,variant_group]
			b = subdf[subdf$Relationship_Label == group2,variant_group]
			
			num_samples = length(a) + length(b)
			
			res = t.test(a,b, paired=TRUE)
			
			pval  = res$p.value
			tstat = res$statistic
			conf  = res$conf.int[1]
			conf2 = res$conf.int[2]
			
			res2 = cohen.d(a, b, paired=TRUE, na.rm=TRUE)			
			cohen_d    = res2$estimate
			
			summ_prs[nrow(summ_prs) + 1,] = c(pval, tstat, conf, conf2, group1, 'different than', group2, variant_group, cohen_d, num_samples)
			
		}
	}	
}

#Concatenate both summ tables together
summ <- rbind(summ_prs,summ_rv)

#---------------------------------
# Filter statistics and bonferroni correction
#---------------------------------


# keep only a few comparisons
comparisons_keep = c(
	'Proband greater than Carrier Parent',
	'Proband greater than Noncarrier Parent',
	'Proband different than Carrier Parent',
	'Proband different than Noncarrier Parent'
)

# create a dummy column to filter comparisons
summ[,'full_comparison'] = glue('{summ$group1} {summ$Comparison} {summ$group2}')

# filter for comparisions in comparisons_keep
summ = summ[summ$full_comparison %in% comparisons_keep,]

# drop the dummy columns
summ$full_comparison = NULL


# FDR correction
summ[,'FDR'] = p.adjust(summ$`P-value`, method='bonferroni')


# file to write out statistics
outfile = 'statistics/t-tests_16p12_max_samples.tsv'

# save
write.table(summ, file=outfile, sep='\t', row.names=F)


#---------------------------------
# heatmap
#---------------------------------

library(corrplot)

make_fake_pmatrix = function(adj_p_mat, unadj_p_mat) {
	# create a fake p value matrix
	# where fdr corrected p < 0.05 gets **
	# and uncorrected p < 0.05 gets *
	fake_p_matrix = adj_p_mat
	for (i in 1:dim(fake_p_matrix)[1]) {
		for (j in 1:dim(fake_p_matrix)[2]) {
			uncorrected_p = unadj_p_mat[i,j]
			corrected_p   = adj_p_mat[i,j]
			if (corrected_p < 0.05) {
				fake_p_matrix[i,j] = 0.002
			} else if (uncorrected_p < 0.05) {
				fake_p_matrix[i,j] = 0.02
			} else {
				fake_p_matrix[i,j] = 0.5
			}

		}
	}
	
	return (fake_p_matrix)
}




df = reshape(summ, idvar = c('Variant Group'), timevar = "group2", direction = "wide")
rownames(df) = df$`Variant Group`

# make p matrix and cohen's d matrix
p_mat = df[c('P-value.Carrier Parent', 'P-value.Noncarrier Parent')]
d_mat = df[c('Cohens D.Carrier Parent', 'Cohens D.Noncarrier Parent')]
adj_p_mat = df[c('FDR.Carrier Parent', 'FDR.Noncarrier Parent')]


cols = c('Carrier Parent', 'Noncarrier Parent')
colnames(d_mat) = cols
colnames(p_mat) = cols
colnames(adj_p_mat) = cols

for (col in cols) {
	d_mat[,col] = as.numeric(d_mat[,col])
	p_mat[,col] = as.numeric(p_mat[,col])
	adj_p_mat[,col] = as.numeric(adj_p_mat[,col])
}

d_mat = as.matrix(d_mat)
p_mat = as.matrix(p_mat)
adj_p_mat = as.matrix(adj_p_mat)
fake_p_mat = make_fake_pmatrix(adj_p_mat, p_mat)

z_max = max(abs(d_mat))

title_text = glue('Burden Tests')
filename = 'figures/Heatmap_max_samples.pdf'
pdf(filename)
corrplot(
	d_mat,
	p.mat = fake_p_mat, 
	method = "color",
	insig = "label_sig",
	sig.level = c(.001, .01, .05),
	pch.cex = .75,
	tl.col = "black",  
	pch.col = "white",
	mar=c(0,0,1,0),
	na.label="-",
	title=title_text,
	cl.ratio = 1,
	col = colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", 
	    "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", 
	    "#D6604D", "#B2182B", "#67001F"))(200),
	    is.corr=F,
	 	cl.lim = c(-z_max, z_max)
)
dev.off()

# col = colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", 
	    # "#FDDBC7", "#F7F7F7", "#D1E5F0", "#92C5DE", 
	    # "#4393C3", "#2166AC", "#053061"))(200)


