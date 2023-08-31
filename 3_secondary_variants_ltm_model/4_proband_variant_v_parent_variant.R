library(effsize)	  # cohen's D
library(glue)

df <- read.csv('../../../11_Variant Integration/16p12_cohort_summary_v22.csv')
df <- df[df$Sample!='', ]

df$Relationship_Label <- ''
df[df$Relationship %in% c('MC', 'FC'), 'Relationship_Label'] <- 'Carrier Parent'
df[df$Relationship %in% c('MNC', 'FNC'), 'Relationship_Label'] <- 'Noncarrier Parent'
df[df$Relationship=='P', 'Relationship_Label'] <- 'Proband'
df <- df[df$Relationship_Label!='', ]

# Remove families with 2 noncarrier parents
nc2_fams <- table(df[df$Relationship_Label=='Noncarrier Parent', 'Family'])
nc2_fams <- rownames(nc2_fams[nc2_fams==2])
df <- df[!(df$Family %in% nc2_fams), ]

# Order by family
df = df[order(df$Relationship_Label),]  
df = df[order(df$Family),]  

sample_groups = c('Proband', 'Carrier Parent', 'Noncarrier Parent')
variant_groups = c('Missense_CADD25', 'Missense_CADD25_LOEUF035', 'LOF', 'LOF_LOEUF_035', 'Splice_CADD25', 'Splice_CADD25_LOEUF035', 'genes_del', 'genes_dup', 'dels_loeuf', 'dups_loeuf', 'STRs_exonic', 'STRs_exonic_LOEUF035', 'Rare_Deleterious_SNVs', 'Rare_Deleterious_SNVs_LOEUF', 'enhancer', 'promoter', 'UTR5', 'intelligence_PRS', 'SCZ_PRS', 'educational_attainment_PRS', 'autism_PRS')

stats_df <- data.frame(Variant1=character(0), Variant2=character(0), group1=character(0), group2=character(0), Pvalue=numeric(0), Num_samples=numeric(0), Coeff=numeric(0))
stats_df$Variant1 = as.character(stats_df$Variant1)
stats_df$Variant2 = as.character(stats_df$Variant2)
stats_df$group1 = as.character(stats_df$group1)
stats_df$group2 = as.character(stats_df$group2)

for (group1 in sample_groups) {
	for (group2 in sample_groups) {
		if (group1 == group2) next
		
		# FILTER
		families_of_group1 = df[df$Relationship_Label == group1,]$Family
		families_of_group2 = df[df$Relationship_Label == group2,]$Family
		families = intersect(families_of_group1, families_of_group2)
		subdf = df[df$Family %in% families,]
		subdf = subdf[subdf$Relationship_Label %in% c(group1, group2), ]

		for (variant_col1 in variant_groups){
			for (variant_col2 in variant_groups){
				
				tmpdf = subdf[c('Family', 'Relationship_Label', variant_col1, variant_col2)]
				tmpdf = na.omit(tmpdf)
				
				# FILTER
				families_of_group1 = tmpdf[tmpdf$Relationship_Label == group1,]$Family
				families_of_group2 = tmpdf[tmpdf$Relationship_Label == group2,]$Family
				families = intersect(families_of_group1, families_of_group2)
				tmpdf = tmpdf[tmpdf$Family %in% families,]

				num_samples = dim(tmpdf)[1]

				a = tmpdf[tmpdf$Relationship_Label == group1, variant_col1]
				b = tmpdf[tmpdf$Relationship_Label == group2,  variant_col2]
				res = cor.test(a, b)
				coeff = as.numeric(res$estimate)
				p_val = res$p.value
				

				app = c(variant_col1, variant_col2, group1, group2, p_val, num_samples, coeff)
				stats_df[nrow(stats_df)+1,] = app


			}
		}
	}
}



#---------------------------------
# Filter statistics and BH correction
#---------------------------------


# keep only a few comparisons
comparisons_keep = c(
	'Proband Carrier Parent',
	'Proband Noncarrier Parent',
	'Carrier Parent Noncarrier Parent'
)

# create a dummy column to filter comparisons
stats_df[,'full_comparison'] = glue('{stats_df$group1} {stats_df$group2}')
stats_df[,'variant_group1'] = glue('{stats_df$Variant1} {stats_df$group1}')
stats_df[,'variant_group2'] = glue('{stats_df$Variant2} {stats_df$group2}')

# filter for comparisions in comparisons_keep
stats_df = stats_df[stats_df$full_comparison %in% comparisons_keep,]



stats_df[,'FDR'] = NA
for (comparison_keep in comparisons_keep) {
	subdf = stats_df[stats_df$full_comparison == comparison_keep,]
	
	pvalues = subdf$Pvalue
	adj_pvalues = p.adjust(pvalues, method='BH')
	subdf[,'FDR'] = adj_pvalues
	
	for (index in rownames(subdf)){
		fdr = subdf[index, 'FDR']
		stats_df[index, 'FDR'] = fdr
	}
}








print(head(stats_df))
print(tail(stats_df))


write.csv(stats_df, 'statistics/variant_v_variant.csv', row.names=F)























