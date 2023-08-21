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

filename = dropbox_location + '/16p12.2 project/Human patients project/WGS paper/11_Variant Integration/16p12_cohort_summary_v18.xlsx'

df = pd.read_excel(filename)

# FILTER 1
df['Relationship_Label'] = df['Relationship'].apply(get_sample_type)
df = df[df.Relationship_Label != 'Other']

# FILTER 2
df = df[~df.Rare_Deleterious_SNVs.isna()]

# FILTER 3
# remove samples if two non carrier parents
families_of_non_carrier_parents = df[df.Relationship_Label == 'Noncarrier Parent'].Family
families_with_two_noncarrier_parents = families_of_non_carrier_parents[families_of_non_carrier_parents.duplicated()].to_list()
df = df[~((df.Family.isin(families_with_two_noncarrier_parents)) & (df.Relationship_Label == 'Noncarrier Parent'))]






def boxplot(df, stats_df, variant_group, relationship):
	order = ['Proband', relationship]
	pairs = [('Proband', relationship)]
	line_offset_to_box = 0.2


	pvalue1 = stats_df[stats_df.group2 == relationship].iloc[0]['P-value']
	
	pvalues=[pvalue1]
	annotations = ['p={:.3}'.format(s) for s in pvalues]
	
	fig = plt.figure(figsize=(4,4))
	
	g = sns.boxplot(data=df, x='Relationship_Label', y=variant_group, whis=[0, 100])
	
	g.set_title(variant_group)
	g.set_xlabel(None)
	
	sns.stripplot(x="Relationship_Label", y=variant_group, data=df,
              size=4, marker="$\circ$", s=5, color='black')
	
	add_stat_annotation(
		g, data=df, x='Relationship_Label', y=variant_group, order=order,
		box_pairs=pairs,
		line_offset_to_box=line_offset_to_box,
		perform_stat_test=False,
		text_format='simple',
		pvalues=pvalues,
		text_annot_custom=annotations
      )
	
	g.set_xlabel(None)
	
	pdf.savefig(fig, bbox_inches='tight')
	plt.close()
	
	


#---------------------------------
# Coding and noncoding variants
#---------------------------------


stats_df = pd.read_csv('statistics/t-tests_16p12_max_samples.tsv', sep='\t')

variant_groups = ['Missense_CADD25', 'Missense_CADD25_LOEUF035', 'LOF', 'LOF_LOEUF_035', 'Splice_CADD25', 'Splice_CADD25_LOEUF035', 'genes_del', 'genes_dup', 'dels_loeuf', 'dups_loeuf', 'STRs_exonic', 'STRs_exonic_LOEUF035', 'Rare_Deleterious_SNVs', 'Rare_Deleterious_SNVs_LOEUF', 'enhancer','promoter','UTR5']

relationships = ['Carrier Parent', 'Noncarrier Parent']


outfile = 'figures/WGS Boxplots.pdf'
pdf = PdfPages(outfile)
for variant_group in variant_groups:
	for relationship in relationships:
		print(variant_group)
		subdf = df[(~df[variant_group].isna()) & (df.Relationship_Label.isin(['Proband', relationship]))]
		sub_stats_df = stats_df[(stats_df['Variant Group'] == variant_group) & (stats_df['group2'] == relationship)]
		
		boxplot(subdf, sub_stats_df, variant_group, relationship)
	


pdf.close()



#---------------------------------
# PRS
#---------------------------------

filename = dropbox_location + '/16p12.2 project/Human patients project/WGS paper/11_Variant Integration/16p12_cohort_summary_v17.xlsx'

df = pd.read_excel(filename)

# FILTER 1
df['Relationship_Label'] = df['Relationship'].apply(get_sample_type)
df = df[df.Relationship_Label != 'Other']

# FILTER 2
df = df[~df.SCZ_PRS.isna()]

# FILTER 3
# remove samples if two non carrier parents
families_of_non_carrier_parents = df[df.Relationship_Label == 'Noncarrier Parent'].Family
families_with_two_noncarrier_parents = families_of_non_carrier_parents[families_of_non_carrier_parents.duplicated()].to_list()
df = df[~((df.Family.isin(families_with_two_noncarrier_parents)) & (df.Relationship_Label == 'Noncarrier Parent'))]




stats_df = pd.read_csv('statistics/t-tests_16p12_max_samples.tsv', sep='\t')

variant_groups = ['intelligence_PRS','SCZ_PRS','educational_attainment_PRS','autism_PRS']


outfile = 'figures/PRS boxplots.pdf'
pdf = PdfPages(outfile)
for variant_group in variant_groups:
	for relationship in relationships:
		print(variant_group)
		subdf = df[(~df[variant_group].isna()) & (df.Relationship_Label.isin(['Proband', relationship]))]
		sub_stats_df = stats_df[(stats_df['Variant Group'] == variant_group) & (stats_df['group2'] == relationship)]
		
		boxplot(subdf, sub_stats_df, variant_group, relationship)
	



pdf.close()
