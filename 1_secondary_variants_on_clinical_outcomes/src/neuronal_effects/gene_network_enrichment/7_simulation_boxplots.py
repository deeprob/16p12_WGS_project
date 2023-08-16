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
from matplotlib.backends.backend_pdf import PdfPages


phenotypes = ['Child_ID_DD', 'Child_behav', 'Child_psych', 'Child_nervous_system', 'Child_congenital', 'Child_craniofacial']
degree_bins = ['0-18', '19-95', '96-329', '330-Max']

df = pd.read_csv('statistics/simulations_by_phenotype.csv')


pdf = PdfPages('figures/simulation_boxplots_by_phenotype.pdf')



# for variant_class in variant_classes:
for phenotype in phenotypes:
	subdf = df[(df['Phenotype'] == phenotype)]
	subdf = subdf.melt(id_vars=['Variant_class', 'Phenotype', 'Simulation'])
	
	
	true_network = subdf[subdf['Simulation'] == 'True']
	true_network = true_network.set_index('variable')
	subdf = subdf[subdf['Simulation'] != 'True']
		
	fig = plt.figure(figsize=(3.5,4))	
	
	
	g = sns.violinplot(data=subdf, x='variable', y='value', cut=0)
	
	g.set_xlabel('Degree bin')
	g.set_ylabel('Number of gene in bin')
	g.set_title(phenotype)
	
	offset = .25
		
	value = true_network.loc['0-18']['value']
	center = 0
	plt.plot([center-offset, center+offset], [value,value], ls='--', color='black')
	
	value = true_network.loc['19-95']['value']
	center = 1
	plt.plot([center-offset, center+offset], [value,value], ls='--', color='black')
	
	value = true_network.loc['96-329']['value']
	center = 2
	plt.plot([center-offset, center+offset], [value,value], ls='--', color='black')
	
	value = true_network.loc['330-Max']['value']
	center = 3
	plt.plot([center-offset, center+offset], [value,value], ls='--', color='black')
	
	
	pdf.savefig(fig, bbox_inches='tight')
	plt.close()



pdf.close()



# ================
# all cohorts
# ================





df = pd.read_csv('statistics/simulations_all_cohorts.csv')



pdf = PdfPages('figures/simulation_boxplots_all_cohorts.pdf')



# for variant_class in variant_classes:
for cohort in df.Cohort.unique():
	for first_hit in df.First_hit.unique():
		subdf = df[(df['Cohort'] == cohort) & (df['First_hit'] == first_hit)]
		
		print(cohort, first_hit)
		if subdf.shape[0] == 0:
			continue
				

		subdf = subdf.melt(id_vars=['Variant_class', 'Cohort', 'First_hit', 'Simulation'])
	
	
		true_network = subdf[subdf['Simulation'] == 'True']
		true_network = true_network.set_index('variable')
		subdf = subdf[subdf['Simulation'] != 'True']
		
		fig = plt.figure(figsize=(3.5,4))	
		
		
		g = sns.violinplot(data=subdf, x='variable', y='value', cut=0)
		
		g.set_xlabel('Degree bin')
		g.set_ylabel('Number of gene in bin')
		g.set_title(cohort + ' ' + first_hit)
		
		offset = .25
			
		value = true_network.loc['0-18']['value']
		center = 0
		plt.plot([center-offset, center+offset], [value,value], ls='--', color='black')
		
		value = true_network.loc['19-95']['value']
		center = 1
		plt.plot([center-offset, center+offset], [value,value], ls='--', color='black')
		
		value = true_network.loc['96-329']['value']
		center = 2
		plt.plot([center-offset, center+offset], [value,value], ls='--', color='black')
		
		value = true_network.loc['330-Max']['value']
		center = 3
		plt.plot([center-offset, center+offset], [value,value], ls='--', color='black')
		
		
		pdf.savefig(fig, bbox_inches='tight')
		plt.close()



pdf.close()


