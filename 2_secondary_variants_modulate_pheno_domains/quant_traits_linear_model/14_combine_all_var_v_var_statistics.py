#!/bin/python3



import pandas as pd


df = pd.DataFrame()


# variant vs. variant
filename = 'statistics/variant_v_variant.csv'
app = pd.read_csv(filename)
app.columns = ['Variant class/PRS1', 'Variant class/PRS2', 'Sample_Group', 'Pvalue', 'Sample size', 'Direction', 'use_in_fdr_calculation', 'FDR']
app['Direction metric'] = 'Pearson\'s R'
app['Test'] = 'Pearson\'s R'
app['Other info'] = ''
app['Group'] = '16p12_WGS'
df = df.append(app)

# variant vs. variant by phenotypic domain
filename = 'statistics/variant_v_variant_phenotypic_domains.csv'
app = pd.read_csv(filename)
app.columns = ['Variant class/PRS1', 'Variant class/PRS2', 'Sample_Group', 'Pvalue', 'Sample size', 'Direction', 'use_in_fdr_calculation', 'FDR']
app['Direction metric'] = 'Cohen\'s D'
app['Test'] = 'T-test'
app['Other info'] = ''
app['Group'] = '16p12_WGS'
df = df.append(app)



       
outfile = 'statistics/all_statistics_variants.xlsx'
df.to_excel(outfile, index=False)


