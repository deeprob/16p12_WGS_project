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





			


filename = '../../../11_Variant Integration/16p12_cohort_summary_v21.xlsx'
df = pd.read_excel(filename)


# FILTERs
df = df[df.Relationship == 'P']

df['All_rare_del_var'] = df['Missense_CADD25'] + df['LOF'] + df['Splice_CADD25'] + df[ 'genes_del'] + df[ 'genes_dup'] +	df['STRs_exonic']

filename = 'statistics/correlation_phenotypes.csv'
stats_df = pd.read_csv(filename, sep=',')
stats_df = stats_df.set_index(['Variant', 'Phenotype'], drop=False)

phenotypes = ['bmi_z_score', 'head_circumference_z_score', 'hrs_mat', 'srs']

def linear_reg(x, y):
	x = np.array(x).reshape(-1,1)
	reg = LinearRegression().fit(x, y)
	return reg.coef_[0], reg.intercept_

outfile = 'figures/scatterplots.pdf'
pdf = PdfPages(outfile)

variants = ['LOF_LOEUF_035','dups_loeuf','STRs_exonic','intelligence_PRS','educational_attainment_PRS']
for variant in variants:

	for pheno in phenotypes:
		print(pheno)
		# for prs_score in prs_scores:
		
		cases_df = df[~df[variant].isna()]
		cases_df = cases_df[~cases_df[pheno].isna()]
		
		case_R = stats_df.loc[(variant, pheno)]['Coeff']

		case_p = stats_df.loc[(variant, pheno)]['Pvalue']

		case_fdr = stats_df.loc[(variant, pheno)]['FDR']

		# fit a line to the scatter plot
		case_coeff, case_y = linear_reg(cases_df[variant], cases_df[pheno])
		print(case_y, case_coeff)
		
		
		fig, ax = plt.subplots(nrows=1,ncols=1, figsize=(3,3))
		g1 = sns.scatterplot(data=cases_df, x=variant, y = pheno)


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
		
		








