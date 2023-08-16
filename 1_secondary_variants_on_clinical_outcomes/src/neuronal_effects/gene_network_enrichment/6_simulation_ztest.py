#!/bin/python3




import pandas as pd
import numpy as np
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests


degree_bins = ['0-18', '19-95', '96-329', '330-Max']
variant_class = 'all_variants'

#======================
# by phenotype
#======================

phenotypes = ['Child_ID_DD', 'Child_behav', 'Child_psych', 'Child_nervous_system', 'Child_congenital', 'Child_craniofacial']


df = pd.read_csv('statistics/simulations_by_phenotype.csv')


summ = []
for phenotype in phenotypes:
	subdf = df[(df['Phenotype'] == phenotype)]

	true_network = subdf[subdf['Simulation'] == 'True']
	true_network = true_network.iloc[0]
	subdf = subdf[subdf['Simulation'] != 'True']
	
	for degree_bin in degree_bins:
		true_value = true_network[degree_bin]
		simulation_values = subdf[degree_bin].to_list()
		
		mean_simulations_values = np.mean(simulation_values)
		std_dev_simulations_values = np.std(simulation_values)

		zscore = (true_value - mean_simulations_values) / std_dev_simulations_values
		pval = norm.sf(abs(zscore))
		
		app = [variant_class, phenotype, degree_bin, zscore, pval, mean_simulations_values, std_dev_simulations_values, true_value]
		summ.append(app)


summ_df = pd.DataFrame(summ, columns=['Variant_class', 'Phenotype', 'Degree_bin', 'Z-score', 'P-value', 'Mean_shortest_distance', 'Std_dev_shortest_distance', 'True_shortest_distance'])



summ_df['FDR'] = multipletests(summ_df['P-value'], method='fdr_bh')[1]

summ_df.to_csv('statistics/simulation_by_phenotype_z_scores.csv', index=False)


#======================
# all cohorts
#======================

phenotypes = ['Child_ID_DD', 'Child_behav', 'Child_psych', 'Child_nervous_system', 'Child_congenital', 'Child_craniofacial']


df = pd.read_csv('statistics/simulations_all_cohorts.csv')


summ = []
for cohort in df.Cohort.unique():
	for first_hit in df.First_hit.unique():
		subdf = df[(df['Cohort'] == cohort) & (df['First_hit'] == first_hit)]
		
		print(cohort, first_hit)
		if subdf.shape[0] == 0:
			continue
		
		true_network = subdf[subdf['Simulation'] == 'True']
		true_network = true_network.iloc[0]
		subdf = subdf[subdf['Simulation'] != 'True']
		
		for degree_bin in degree_bins:
			true_value = true_network[degree_bin]
			simulation_values = subdf[degree_bin].to_list()
			
			mean_simulations_values = np.mean(simulation_values)
			std_dev_simulations_values = np.std(simulation_values)

			zscore = (true_value - mean_simulations_values) / std_dev_simulations_values
			pval = norm.sf(abs(zscore))
			
			app = [variant_class, cohort, first_hit, degree_bin, zscore, pval, mean_simulations_values, std_dev_simulations_values, true_value]
			summ.append(app)


summ_df = pd.DataFrame(summ, columns=['Variant_class', 'Cohort', 'First_hit', 'Degree_bin', 'Z-score', 'P-value', 'Mean_shortest_distance', 'Std_dev_shortest_distance', 'True_shortest_distance'])



summ_df['FDR'] = multipletests(summ_df['P-value'], method='fdr_bh')[1]

summ_df.to_csv('statistics/simulation_all_cohorts_z_scores.csv', index=False)















