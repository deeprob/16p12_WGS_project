



library(corrplot)     # correlation plots
library(reshape2)     # to melt the dataframe
library(glue)



df = read.csv('statistics/variant_v_variant_phenotypic_domains.csv', check.names=F)




phenotype_groups = c('Child_ID_DD', 'Child_behav', 'Child_psych', 'Child_nervous_system', 'Child_congenital', 'Child_craniofacial')
variant_groups = c('Missense_CADD25', 'Missense_CADD25_LOEUF035', 'LOF', 'LOF_LOEUF_035', 'Splice_CADD25', 'Splice_CADD25_LOEUF035', 'genes_del', 'genes_dup', 'dels_loeuf', 'dups_loeuf', 'STRs_exonic', 'STRs_exonic_LOEUF035', 'Rare_Deleterious_SNVs', 'Rare_Deleterious_SNVs_LOEUF', 'enhancer', 'promoter', 'UTR5', 'intelligence_PRS', 'SCZ_PRS', 'educational_attainment_PRS', 'autism_PRS')



make_square_matrix = function(df, column) {
	cor_mat = matrix(nrow=length(variant_groups), ncol=length(variant_groups))
	rownames(cor_mat) = variant_groups
	colnames(cor_mat) = variant_groups
	
	for (col1 in variant_groups) {
		for (col2 in variant_groups) {
			
			value = df[(df$`Variant1` == col1) & (df$Variant2 == col2), column]
			cor_mat[col1, col2] = value
		}
	}
	
	return(cor_mat)
}

make_fake_pmatrix = function(adj_p_mat, unadj_p_mat) {
	# create a fake p value matrix
	# where fdr corrected p < 0.05 gets **
	# and uncorrected p < 0.05 gets *
	fake_p_matrix = adj_p_mat
	for (i in 1:dim(fake_p_matrix)[1]) {
		for (j in 1:dim(fake_p_matrix)[2]) {
			uncorrected_p = unadj_p_mat[i,j]
			corrected_p   = adj_p_mat[i,j]
			if (is.na(corrected_p)) {
				fake_p_matrix[i,j] = 0
			} else if (corrected_p < 0.05) {
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

pdf(glue('figures/correlation_variants_by_phenotype.pdf'))


for (phenotype_group in phenotype_groups) {
	print(phenotype_group)

	
	dfplot = df[df$Phenotype_Group == phenotype_group,]
	min_num_samples = min(dfplot$Num_samples)
	max_num_samples = max(dfplot$Num_samples)
	title_text = glue('{phenotype_group} ({min_num_samples}-{max_num_samples} Samples)')
	p_mat = make_square_matrix(dfplot, 'Pvalue')
	cor_mat = make_square_matrix(dfplot, 'Coeff')

	adj_p_mat = make_square_matrix(dfplot, 'FDR')
	fake_p_mat = make_fake_pmatrix(adj_p_mat, p_mat)


	corrplot(
			cor_mat,
			p.mat = fake_p_mat, 
			method = "circle",
		 	insig = "label_sig",
		 	sig.level = c(.001, .01, .05),
		 	pch.cex = .75,
		 	tl.col = "black",  
		 	pch.col = "white",
			mar=c(0,0,1,0),
			na.label="-",
		 	title=title_text,
		 	col = colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))(20),
		 	type='upper',
		)

}







dev.off()
