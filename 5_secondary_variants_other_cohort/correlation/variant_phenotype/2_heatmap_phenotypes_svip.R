



library(corrplot)     # correlation plots
library(reshape2)     # to melt the dataframe
library(glue)



df = read.csv('statistics/variant_v_phenotype.csv', check.names=F)


print(head(df))

groups = c('16p-deletion', '16p-duplication')
variant_cols = c('Missense CADD>25', 'CADD25 LOEUF<0.35', 'LOF CADD25_NA', 'LOF CADD LOEUF<0.35', 'Splice CADD25', 'Splice CADD25 LOEUF<0.35', 'Rare_Deleterious_SNVs', 'Rare_Deleterious_SNVs_LOEUF', 'Genes DEL', 'Genes DUP', 'DELs LOEUF<0.35', 'DUPs LOEUF<0.35', 'intelligence_PRS', 'SCZ_PRS', 'educational_attainment_PRS', 'autism_PRS')
pheno_cols = c('Full-scale IQ','SRS raw','ABCL_CBCL external','ABCL_CBCL internal','BSI','BMI Z-score','Head circumference Z-score')



make_square_matrix = function(df, column) {
	cor_mat = matrix(nrow=length(variant_cols), ncol=length(pheno_cols))
	rownames(cor_mat) = variant_cols
	colnames(cor_mat) = pheno_cols
	
	
	for (col1 in variant_cols) {
		for (col2 in pheno_cols) { 
			value = df[(df$`Variant class/PRS` == col1) & (df$Phenotype == col2), column]
			
			
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
				fake_p_matrix[i,j] = NA
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

pdf(glue('figures/correlation_phenotypes.pdf'))


for (group in groups) {
	print(group)

	
	dfplot = df[df$Group == group,]
	min_num_samples = min(dfplot$`Sample size`)
	max_num_samples = max(dfplot$`Sample size`)
	title_text = glue('{group} ({min_num_samples}-{max_num_samples} Probands)')
	p_mat = make_square_matrix(dfplot, 'P-value')
	cor_mat = make_square_matrix(dfplot, 'Direction')
	adj_p_mat = make_square_matrix(dfplot, 'FDR')
	fake_p_mat = make_fake_pmatrix(adj_p_mat, p_mat)


	corrplot(
			cor_mat,
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
		 	col = colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))(20)
		)

}







dev.off()
