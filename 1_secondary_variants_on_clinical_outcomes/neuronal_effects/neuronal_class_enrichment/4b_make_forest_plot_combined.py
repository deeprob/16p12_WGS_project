#!/bin/python3


import pandas as pd
import sys
import numpy as np


# libraries related to plotting
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
sns.set_style({'font.family':'sans-serif', 'font.sans-serif':'Arial'})
from matplotlib.backends.backend_pdf import PdfPages




super_df = pd.read_csv('statistics/fishers_exact_with_R_superclass.csv')
subclass_df = pd.read_csv('statistics/fishers_exact_with_R_subclass.csv')
df = pd.concat([super_df,subclass_df])
#meta_df = pd.read_csv('allen_institute_tables/metadata.csv')
#meta_df = meta_df.set_index('cluster_label', drop=False)




# only keep significant
df = df[df.FDR <= 0.05]

# add class and subclass label
#df['cluster_label'] = df['cell_type'].apply(lambda s: s.split('|')[0])
#cluster_label_dict = meta_df[['cluster_label', 'class_label']].drop_duplicates()
#df['class_label'] = df['cluster_label'].map(cluster_label_dict['class_label'])
#cluster_label_dict = meta_df[['cluster_label', 'subclass_label']].drop_duplicates()
#df['subclass_label'] = df['cluster_label'].map(cluster_label_dict['subclass_label'])


# will use subclass label to plot the statistics
# drop the duplicate subclass labels keeping the highest oddsratio
subdf = df.sort_values('oddsratio', ascending=False)
#subdf = df.drop_duplicates(['class_label', 'subclass_label'])
#subdf = subdf.sort_values(['class_label', 'subclass_label'])
subdf['label'] = subdf['cell_type']+"_"+subdf['variant_class']


# reverse order and reset index
#subdf = subdf.iloc[::-1]
subdf = subdf.reset_index()



pdf = PdfPages('figures/forest_plot_combined.pdf')
fig = plt.figure(figsize=(2.5,4))

# draw points
x = subdf['oddsratio']
y = subdf['label']
# plt.scatter(x,y)
plt.scatter(x,y, color='tab:green')

# draw horizontal lines
for i in subdf.index:
	lower = subdf.loc[i, 'conf_int_lower']
	upper = subdf.loc[i, 'conf_int_upper']
	
	x = [lower, upper]
	y = [i, i]
	# plt.plot(x,y, color='#1f77b4')
	plt.plot(x, y, color='tab:green')


# draw vertical line
plt.axvline(x=1.0, color='black', ls='--')

# axes parameters
plt.xlim([0,3])
ax = plt.gca()
ax.xaxis.set_minor_locator(MultipleLocator(0.1))


pdf.savefig(fig, bbox_inches='tight')

pdf.close()

exit()

_________________________



