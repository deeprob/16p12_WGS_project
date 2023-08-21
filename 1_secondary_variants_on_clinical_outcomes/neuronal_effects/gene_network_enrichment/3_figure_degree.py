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
from statannot import add_stat_annotation
from matplotlib.backends.backend_pdf import PdfPages





def plot(df, labels_order):
	fig = plt.figure(figsize=(3.5,4))
	g = sns.boxplot(data=df, x='Quartile', y='LOEUF', order=labels_order, whis=[0,100])
	g.set_title('Brain-Specific Network')
	g.set_xlabel('Node degree bin')
	# g.set_ylim([0,2])

	pdf.savefig(fig, bbox_inches='tight')
	
	plt.close()

def plot2(df):
	fig = plt.figure()
	g = sns.scatterplot(data=df, x='Degree', y='LOEUF')
	g.set_title('Brain-Specific Network')

	pdf.savefig(fig, bbox_inches='tight')
	
	plt.close()

def plot3(df):
	fig = plt.figure()
	g = sns.histplot(data=df, x='LOEUF')
	g.set_title('Brain-Specific Network')
	pdf.savefig(fig, bbox_inches='tight')
	plt.close()

def plot4(df):
	fig = plt.figure()
	g = sns.histplot(data=df, x='Degree')
	g.set_title('Brain-Specific Network')
	pdf.savefig(fig, bbox_inches='tight')
	plt.close()




pdf = PdfPages('figures/Degree_Loeuf.pdf')
df = pd.read_csv('statistics/network_degree_loeuf.csv')
df = df[~df.LOEUF.isna()]
df.LOEUF = df.LOEUF.astype(float)

labels_order = ['0-18','19-95','96-329','330-Max']
plot(df, labels_order)
plot2(df)
plot3(df)
plot4(df)

pdf.close()
	



