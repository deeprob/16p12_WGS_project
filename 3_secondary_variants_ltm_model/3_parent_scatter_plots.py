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




dropbox_location = 'C:/Users/corny/Dropbox/'
			
def get_sample_type(s):
	if s == 'P':
		return 'Proband'
	elif s == 'MC':
		return 'Carrier_parent'
	elif s == 'FC':
		return 'Carrier_parent'
	elif s == 'FNC':
		return 'Noncarrier_parent'
	elif s == 'MNC':
		return 'Noncarrier_parent'
	else:
		return 'Other'

filename = dropbox_location+'16p12.2 project/Human patients project/WGS paper/11_Variant Integration/16p12_cohort_summary_v22.csv'
df = pd.read_csv(filename)


df['Relationship_Label'] = df['Relationship'].apply(get_sample_type)
df = df[df.Relationship_Label.isin(['Proband', 'Carrier_parent', 'Noncarrier_parent'])]
df = df[~df.SCZ_PRS.isnull()]
df = df[~df['Missense_CADD25'].isnull()]

df['All_rare_del_var'] = df['Missense_CADD25'] + df['LOF'] + df['Splice_CADD25'] + df[ 'genes_del'] + df[ 'genes_dup'] + df['STRs_exonic']

filename = 'statistics/compared_to_parents_correlations.tsv'
stats_df = pd.read_csv(filename, sep='\t')
stats_df = stats_df.set_index(['Role', 'PRS'], drop=False)

prs_scores= ['SCZ_PRS', 'intelligence_PRS', 'educational_attainment_PRS', 'autism_PRS']




def linear_reg(x, y):
	x = np.array(x).reshape(-1,1)
	reg = LinearRegression().fit(x, y)
	return reg.coef_[0], reg.intercept_

outfile = 'figures/scatterplots_parents_share_axis.pdf'
pdf = PdfPages(outfile)

for prs_score in prs_scores:
	proband_df = df[df['Relationship_Label'] == 'Proband']
	carrier_parent_df = df[df['Relationship_Label'] == 'Carrier_parent']
	noncarrier_parent_df = df[df['Relationship_Label'] == 'Noncarrier_parent']


	proband_R = stats_df.loc[('Proband', prs_score)]['R']
	carrier_parent_R = stats_df.loc[('Carrier_parent', prs_score)]['R']
	noncarrier_parent_R = stats_df.loc[('Noncarrier_parent', prs_score)]['R']

	proband_p = stats_df.loc[('Proband', prs_score)]['pvalue']
	carrier_parent_p = stats_df.loc[('Carrier_parent', prs_score)]['pvalue']
	noncarrier_parent_p = stats_df.loc[('Noncarrier_parent', prs_score)]['pvalue']

	proband_fdr = stats_df.loc[('Proband', prs_score)]['FDR']
	carrier_parent_fdr = stats_df.loc[('Carrier_parent', prs_score)]['FDR']
	noncarrier_parent_fdr = stats_df.loc[('Noncarrier_parent', prs_score)]['FDR']


	proband_coeff, proband_y = linear_reg(proband_df[prs_score], proband_df['All_rare_del_var'])
	carrier_parent_coeff, carrier_parent_y = linear_reg(carrier_parent_df[prs_score], carrier_parent_df['All_rare_del_var'])
	noncarrier_parent_coeff, noncarrier_parent_y = linear_reg(noncarrier_parent_df[prs_score], noncarrier_parent_df['All_rare_del_var'])


	fig, ax = plt.subplots(nrows=1,ncols=3, figsize=(8,3), sharex=True, sharey=True)
	g1 = sns.scatterplot(data=proband_df, x=prs_score, y = 'All_rare_del_var', ax=ax[0])
	g2 = sns.scatterplot(data=carrier_parent_df, x=prs_score, y = 'All_rare_del_var', ax=ax[1])
	g3 = sns.scatterplot(data=noncarrier_parent_df, x=prs_score, y = 'All_rare_del_var', ax=ax[2])


	text = 'R = {:.2f}\np = {:.2f}\nFDR = {:.2f}'.format(proband_R, proband_p, proband_fdr)
	g1.annotate(text, xy=(0.58, 0.75), xycoords='axes fraction')
	text = 'R = {:.2f}\np = {:.2f}\nFDR = {:.2f}'.format(carrier_parent_R, carrier_parent_p, carrier_parent_fdr)
	g2.annotate(text, xy=(0.58, 0.75), xycoords='axes fraction')
	text = 'R = {:.2f}\np = {:.2f}\nFDR = {:.2f}'.format(noncarrier_parent_R, noncarrier_parent_p, noncarrier_parent_fdr)
	g3.annotate(text, xy=(0.58, 0.75), xycoords='axes fraction')


	x_vals = np.array(g1.get_xlim())
	y_vals = proband_y + proband_coeff * x_vals
	sns.lineplot(x=x_vals, y=y_vals, ax=ax[0])
	x_vals = np.array(g2.get_xlim())
	y_vals = carrier_parent_y + carrier_parent_coeff * x_vals
	sns.lineplot(x=x_vals, y=y_vals, ax=ax[1])
	x_vals = np.array(g3.get_xlim())
	y_vals = noncarrier_parent_y + noncarrier_parent_coeff * x_vals
	sns.lineplot(x=x_vals, y=y_vals, ax=ax[2])
	# print(x_vals, y_vals)
	# x_vals = np.array(g2.get_xlim())
	# y_vals = control_y + control_coeff * x_vals
	# sns.lineplot(x_vals, y_vals, ax=ax[1])



	g1.set_title('Probands')
	g2.set_title('Carrier Parents')
	g3.set_title('Noncarrier Parents')

	plt.tight_layout()
	pdf.savefig(fig, bbox_inches='tight')
	plt.close()

		

pdf.close()
		
		








