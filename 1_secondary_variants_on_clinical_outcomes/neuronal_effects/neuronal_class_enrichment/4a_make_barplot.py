#!/bin/python3




import pandas as pd
import subprocess

# libraries related to plotting
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
sns.set_style({'font.family':'sans-serif', 'font.sans-serif':'Arial'})
from matplotlib.backends.backend_pdf import PdfPages


subprocess.run('mkdir figures', shell=True)


df = pd.read_csv('statistics/fishers_exact.csv')
meta_df = pd.read_csv('allen_institute_tables/metadata.csv')
meta_df = meta_df.set_index('cluster_label', drop=False)




# only keep significant
df = df[df.FDR <= 0.05]

# add class and subclass label
df['cluster_label'] = df['cell_type'].apply(lambda s: s.split('|')[0])
cluster_label_dict = meta_df[['cluster_label', 'class_label']].drop_duplicates()
df['class_label'] = df['cluster_label'].map(cluster_label_dict['class_label'])
cluster_label_dict = meta_df[['cluster_label', 'subclass_label']].drop_duplicates()
df['subclass_label'] = df['cluster_label'].map(cluster_label_dict['subclass_label'])


# will use subclass label to plot the statistics
# drop the duplicate subclass labels keeping the highest oddsratio
df = df.sort_values('oddsratio', ascending=False)
subdf = df.drop_duplicates(['class_label', 'subclass_label'])
subdf = subdf.sort_values(['class_label', 'subclass_label'])
subdf['label'] = subdf['class_label'] + '|' + subdf['subclass_label']


print(df)
print(subdf)


pdf = PdfPages('figures/barplot.pdf')
fig = plt.figure(figsize=(2.5,4))

g = sns.barplot(data=subdf, x='oddsratio', y='subclass_label', color='tab:blue')
plt.axvline(x=1.0, color='black', ls='--')

g.set_ylabel(None)
g.set_xlim([0,3])
ax = plt.gca()
ax.xaxis.set_minor_locator(MultipleLocator(0.1))


pdf.savefig(fig, bbox_inches='tight')

pdf.close()

































