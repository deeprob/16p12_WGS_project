#!/bin/R



library(corrplot)
library(glue)

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

make_square_matrix = function(df, group1, group2, column) {
	cor_mat = matrix(nrow=length(variant_groups), ncol=length(variant_groups))
	rownames(cor_mat) = paste(variant_groups, group1, sep=' ')
	colnames(cor_mat) = paste(variant_groups, group2, sep=' ')
	
	for (col1 in variant_groups) {
		for (col2 in variant_groups) {
			
			value = df[(df$variant_group1 == glue('{col1} {group1}')) & (df$variant_group2 == glue('{col2} {group2}')), column]
			cor_mat[glue('{col1} {group1}'), glue('{col2} {group2}')] = value
		}
	}
	
	return(cor_mat)
}




df = read.csv('statistics/variant_v_variant.csv', check.names=F)


comparisons_keep = list(
	c('Proband', 'Carrier Parent'),
	c('Proband', 'Noncarrier Parent'),
	c('Carrier Parent', 'Noncarrier Parent')
)

for (comparison in comparisons_keep) {
  variant_groups = c('Missense_CADD25', 'Missense_CADD25_LOEUF035', 'LOF', 'LOF_LOEUF_035', 'Splice_CADD25', 'Splice_CADD25_LOEUF035', 'genes_del', 'genes_dup', 'dels_loeuf', 'dups_loeuf', 'STRs_exonic', 'STRs_exonic_LOEUF035', 'Rare_Deleterious_SNVs', 'Rare_Deleterious_SNVs_LOEUF', 'enhancer', 'promoter', 'UTR5', 'intelligence_PRS', 'SCZ_PRS', 'educational_attainment_PRS', 'autism_PRS')
	group1 = comparison[1]
	group2 = comparison[2]
	
	subdf = df[df$full_comparison == glue('{group1} {group2}'), ]


	p_mat = make_square_matrix(subdf, group1, group2, 'Pvalue')
	cor_mat = make_square_matrix(subdf, group1, group2, 'Coeff')

	adj_p_mat = make_square_matrix(subdf, group1, group2, 'FDR')
	fake_p_mat = make_fake_pmatrix(adj_p_mat, p_mat)

	title_text = glue('Burden Tests {group1} vs. {group2}')
	filename = glue('figures/Heatmap_{group1}_{group2}_circle.pdf')
	pdf(filename)
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
		 	col = colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))(200)
		)
	
	# Remove LOUEF variants
	variant_groups = c('Missense_CADD25', 'LOF', 'Splice_CADD25', 'genes_del', 'genes_dup', 'STRs_exonic', 'Rare_Deleterious_SNVs', 'Rare_Deleterious_SNVs_LOEUF', 'enhancer', 'promoter', 'UTR5', 'intelligence_PRS', 'SCZ_PRS', 'educational_attainment_PRS', 'autism_PRS')
	
	p_mat = make_square_matrix(subdf, group1, group2, 'Pvalue')
	cor_mat = make_square_matrix(subdf, group1, group2, 'Coeff')
	  
	adj_p_mat = make_square_matrix(subdf, group1, group2, 'FDR')
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
	 col = colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))(200)
	)
	dev.off()
}


