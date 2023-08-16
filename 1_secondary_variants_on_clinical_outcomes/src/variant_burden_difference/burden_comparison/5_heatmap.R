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





df = read.table('statistics/wilcox_test.tsv', sep='\t', header=T, check.names=F)





df = reshape(df, idvar = c('Variant Group'), timevar = "group2", direction = "wide")
rownames(df) = df$`Variant Group`


cols = c('Proband', 'Carrier Parent')

# make p matrix and cohen's d matrix
p_mat = df[paste('P-value.', cols, sep='')]
d_mat = df[paste('Cohens D.', cols, sep='')]
adj_p_mat = df[paste('FDR.', cols, sep='')]



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
filename = 'figures/Heatmap.pdf'
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
	    "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582", 
	    "#D6604D", "#B2182B", "#67001F"))(20),
	    is.corr=F,
	 	col.lim = c(-z_max, z_max)
)
dev.off()



