#!/bin/python3




import pandas as pd
from sklearn.linear_model import LinearRegression
import numpy as np

# libraries related to plotting
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
sns.set_style({'font.family':'sans-serif', 'font.sans-serif':'Arial'})
from statannot import add_stat_annotation
from matplotlib.backends.backend_pdf import PdfPages




dropbox_location = '~/Dropbox/'
			


filename = dropbox_location+'16p12.2 project/Human patients project/WGS paper/11_Variant Integration/16p12_cohort_summary_v17.xlsx'
df = pd.read_excel(filename)


# FILTERs
df = df[df.Relationship == 'P']
df = df[~df.SCZ_PRS.isna()]
df = df[~df['Missense_CADD25'].isna()]

df['All_rare_del_var'] = df['Missense_CADD25'] + df['LOF'] + df['Splice_CADD25'] + df[ 'genes_del'] + df[ 'genes_dup'] +	df['STRs_exonic']

filename = 'statistics/correlations_rare_variant_prs.tsv'
stats_df = pd.read_csv(filename, sep='\t')
stats_df = stats_df.set_index(['Phenotype', 'PRS'], drop=False)

phenotypes = ['Child_ID_DD', 'Child_behav', 'Child_psych', 'Child_nervous_system', 'Child_congenital', 'Child_craniofacial', 'De_vrie']
prs_scores= ['SCZ_PRS', 'intelligence_PRS', 'educational_attainment_PRS', 'autism_PRS']


phenotype_splits = {
	'Child_ID_DD':2,
	'Child_behav':1,
	'Child_psych':0,
	'Child_nervous_system':0,
	'Child_congenital':1,
	'Child_craniofacial':2,
	'De_vrie':7
}


def linear_reg(x, y):
	x = np.array(x).reshape(-1,1)
	reg = LinearRegression().fit(x, y)
	return reg.coef_[0], reg.intercept_

outfile = 'figures/scatterplots.pdf'
pdf = PdfPages(outfile)

prs_score = 'SCZ_PRS'
for prs_score in prs_scores:
	for pheno in phenotypes:
		print(pheno)
		# for prs_score in prs_scores:
		split_value = phenotype_splits[pheno]
		
		cases_df = df[df[pheno] > split_value]
		
		case_R = stats_df.loc[(pheno, prs_score)]['R']

		case_p = stats_df.loc[(pheno, prs_score)]['pvalue']

		case_fdr = stats_df.loc[(pheno, prs_score)]['FDR']

		# fit a line to the scatter plot
		case_coeff, case_y = linear_reg(cases_df[prs_score], cases_df['All_rare_del_var'])
		print(case_y, case_coeff)
		
		
		fig, ax = plt.subplots(nrows=1,ncols=1, figsize=(3,3))
		g1 = sns.scatterplot(data=cases_df, x=prs_score, y = 'All_rare_del_var')


		text = 'R = {:.2f}\np = {:.2f}\nFDR = {:.2f}'.format(case_R, case_p, case_fdr)
		g1.annotate(text, xy=(0.58, 0.75), xycoords='axes fraction')



		# draw a line that was fit with linear regression
		x_vals = np.array(g1.get_xlim())
		y_vals = case_y + case_coeff * x_vals
		sns.lineplot(x_vals, y_vals)
		print(x_vals, y_vals)




		g1.set_title('Cases')
		fig.suptitle(pheno)

		plt.tight_layout()
		pdf.savefig(fig, bbox_inches='tight')
		plt.close()

		

pdf.close()
		
		








