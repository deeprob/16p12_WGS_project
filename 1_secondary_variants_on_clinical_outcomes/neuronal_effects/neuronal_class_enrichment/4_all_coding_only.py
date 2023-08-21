#!/bin/python3
import pandas as pd
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

# only keep "All coding SNVs"
df=df[df.variant_class=="rare_deleterious_snvs"]
df=df[df.FDR<=0.05]
print(df)

# will use subclass label to plot the statistics
# drop the duplicate subclass labels keeping the highest oddsratio
subdf = df.sort_values('oddsratio', ascending=False)
subdf['label'] = subdf['cell_type']+"_"+subdf['variant_class']

# reset index
subdf = subdf.reset_index()

pdf = PdfPages('figures/4_forest_snvs.pdf')
fig = plt.figure(figsize=(2.5,4))

# draw points
x = subdf['oddsratio']
y = subdf['cell_type']
plt.scatter(x,y, color='tab:green')

# draw horizontal lines
for i in subdf.index:
	lower = subdf.loc[i, 'conf_int_lower']
	upper = subdf.loc[i, 'conf_int_upper']
	
	x = [lower, upper]
	y = [i, i]
	plt.plot(x, y, color='tab:green')

# draw vertical line
plt.axvline(x=1.0, color='black', ls='--')

# axes parameters
plt.xlim([0,2.15])
ax = plt.gca()
ax.xaxis.set_minor_locator(MultipleLocator(0.1))

pdf.savefig(fig, bbox_inches='tight')
pdf.close()