#!/bin/python
import pandas as pd
import numpy as np

# Calculate scores for family history
# Scores will be calculated from parent phenotypes using the adult phenotypic domains
# We will use each of the following methods:
# 1. Max score - the maximum score of both parents for each domain
# 2. Sum score - the sum of the parent scores for each domain
# 3. Mother score - mother's score
# 4. Father score - father's score
# 5. Carrier parent score - carrier parent's score
# 6. NonCarrier parent score - noncarrier parent's score
# Methods 1-4 will only be done on complete trios
# Methods 5 and 6 will only be done on complete trios with an inherited deletion

# Get families
df = pd.read_excel('../11_Variant Integration/16p12_cohort_summary_v19.xlsx')
# Only care about parents with phenotype data
df = df[(df.Relationship.isin(['MC', 'MNC', 'FC', 'FNC', 'M', 'F'])) & (df['Phenotypic Domains']=='X')]
# Restrict to complete families (families with 2 parents)
fam_counts = df.Family.value_counts()
fam_counts = fam_counts[fam_counts == 2]
df = df[df.Family.isin(fam_counts.index.to_list())]

# Reformat data for output
out = pd.DataFrame({'Family':fam_counts.index})
out['Mother'] = out.Family.map(dict(zip(df[df.Relationship.isin(['M', 'MC', 'MNC'])]['Family'].to_list(), df[df.Relationship.isin(['M', 'MC', 'MNC'])]['Sample'].to_list())))
out['Father'] = out.Family.map(dict(zip(df[df.Relationship.isin(['F', 'FC', 'FNC'])]['Family'].to_list(), df[df.Relationship.isin(['F', 'FC', 'FNC'])]['Sample'].to_list())))
out['Carrier'] = out.Family.map(dict(zip(df[df.Relationship.isin(['MC', 'FC'])]['Family'].to_list(), df[df.Relationship.isin(['MC', 'FC'])]['Sample'].to_list())))
out['NonCarrier'] = out.Family.map(dict(zip(df[df.Relationship.isin(['MNC', 'FNC'])]['Family'].to_list(), df[df.Relationship.isin(['MNC', 'FNC'])]['Sample'].to_list())))
out.loc[out.Carrier.isnull(), 'NonCarrier'] = np.nan

# Calculate domain scores for each family for each method
for domain in ['Parent_cognitive', 'Parent_SCZ', 'Parent_depression', 'Parent_addiction']:
    out['max'+domain] = out.Family.apply(lambda fam: max(df[df.Family==fam][domain].to_numpy()))
    out['sum'+domain] = out.Family.apply(lambda fam: sum(df[df.Family==fam][domain].to_numpy()))
    out['mother'+domain] = out.Mother.apply(lambda m: df[df.Sample==m][domain].to_list()[0])
    out['father'+domain] = out.Father.apply(lambda f: df[df.Sample==f][domain].to_list()[0])
    
    out['Carrier'+domain] = np.nan
    out.loc[~out.Carrier.isnull(), 'Carrier'+domain] = out.loc[~out.Carrier.isnull(), 'Carrier'].apply(lambda c: df[df.Sample==c][domain].to_list()[0])
    out['Noncarrier'+domain] = np.nan
    out.loc[~out.NonCarrier.isnull(), 'Noncarrier'+domain] = out.loc[~out.NonCarrier.isnull(), 'NonCarrier'].apply(lambda nc: df[df.Sample==nc][domain].to_list()[0])

# Add probands
df = pd.read_excel('../11_Variant Integration/16p12_cohort_summary_v19.xlsx')
out['Proband'] = out.Family.map(dict(zip(df[df.Relationship=='P']['Family'].to_list(), df[df.Relationship=='P']['Sample'].to_list())))

# Reorganize columns
cols = ['Proband', 'Family', 'Mother', 'Father', 'Carrier', 'NonCarrier']
for method in ['max', 'sum', 'mother', 'father', 'Carrier', 'Noncarrier']:
    for domain in ['Parent_cognitive', 'Parent_SCZ', 'Parent_depression', 'Parent_addiction']:
        cols.append(method+domain)
out = out[cols]

# Save to file
out.to_csv('Analysis_files/1_domain_scores.csv', index = False)