#!/bin/python3




import pandas as pd

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


def get_sample_type(s):
	if (s == 'P'):
		return('Proband')
	elif (s == 'MC'):
		return('Carrier Parent')
	elif (s == 'FC'):
		return('Carrier Parent')
	elif (s == 'FNC'):
		return('Noncarrier Parent')
	elif (s == 'MNC'):
		return('Noncarrier Parent')
	else:
		return('Other')



dropbox_location = '~/Dropbox/'

filename = dropbox_location + '/16p12.2 project/Human patients project/WGS paper/11_Variant Integration/16p12_cohort_summary_v17.xlsx'

df = pd.read_excel(filename)

# FILTER 1
df['Relationship_Label'] = df['Relationship'].apply(get_sample_type)
df = df[df.Relationship_Label.isin(['Proband', 'Carrier Parent'])]


filename = f'{dropbox_location}/16p12.2 project/Human patients project/WGS paper/11_Variant Integration/Estonian_Summary.xlsx'
estonia_df = pd.read_excel(filename)
estonia_df = estonia_df[estonia_df.Relationship == 'P']
estonia_df['Relationship_Label'] = 'Estonia'
df = df.append(estonia_df)


def boxplot(df, variant_group):
	
	fig = plt.figure(figsize=(4,4))
	
	g = sns.boxplot(data=df, x='Relationship_Label', y=variant_group, whis=[0, 100])
	
	g.set_title(variant_group)
	g.set_xlabel(None)
	
	sns.stripplot(x="Relationship_Label", y=variant_group, data=df,
              size=4, marker="$\circ$", s=5, color='black')
	

	g.set_xlabel(None)
	
	pdf.savefig(fig, bbox_inches='tight')
	plt.close()
	
	



variant_groups = ['Missense_CADD25', 'Missense_CADD25_LOEUF035', 'LOF', 'LOF_LOEUF_035', 'Splice_CADD25', 'Splice_CADD25_LOEUF035', 'genes_del', 'genes_dup', 'dels_loeuf', 'dups_loeuf', 'STRs_exonic', 'STRs_exonic_LOEUF035', 'Rare_Deleterious_SNVs', 'Rare_Deleterious_SNVs_LOEUF', 'enhancer','promoter','UTR5','intelligence_PRS','SCZ_PRS','educational_attainment_PRS','autism_PRS']



outfile = 'figures/Boxplots.pdf'
pdf = PdfPages(outfile)
for variant_group in variant_groups:
	print(variant_group)
	subdf = df[~df[variant_group].isna()]
	
	boxplot(subdf, variant_group)
	

pdf.close()


