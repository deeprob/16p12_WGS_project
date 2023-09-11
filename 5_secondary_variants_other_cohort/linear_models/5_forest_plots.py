#!/bin/python3


import pandas as pd
import sys

# libraries related to plotting
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
sns.set_style({'font.family':'sans-serif', 'font.sans-serif':'Arial'})

# name of input file and name of output
infile   = sys.argv[1]
outfile  = sys.argv[2]
title    = sys.argv[3]
width    = float(sys.argv[4])
height   = float(sys.argv[5])

print(infile)
print(outfile)
print(title)

# read in statistics
df = pd.read_csv(infile, sep='\t')

df = df[~df.Variable.isin(['(Intercept)', 'Sex'])]
df = df.reset_index(drop=True)


print(df)
# rename columns
# df.columns = ['Label', 'middle', 'lower', 'upper', 'p-value', 'Num_samples']
print(df.columns)

# add number of samples to title
num_samples = df.iloc[0]['Num_samples']
title = '{} {}'.format(num_samples, title)

# plot
plt.figure(figsize=(width, height))


# draw points
x = df['Log Odds Ratio']
y = df['Variable']

plt.scatter(x,y)

# draw horizontal lines
for i in df.index:
	lower = df.loc[i, '2.5% C.I.']
	upper = df.loc[i, '97.5% C.I.']
	
	x = [lower, upper]
	y = [i, i]
	plt.plot(x,y, color='#1f77b4')

# draw vertical line
x = [0 for s in range(0,len(df.index))]
y = list(df.index)

plt.plot(x, y, linewidth=0.5, color='black', ls='--')

plt.gca().xaxis.set_major_locator(MultipleLocator(1))
plt.gca().xaxis.set_minor_locator(MultipleLocator(.25))

plt.xlabel('Beta Coefficient')

plt.title(title)

# plt.tight_layout()
plt.savefig(outfile, bbox_inches="tight")

