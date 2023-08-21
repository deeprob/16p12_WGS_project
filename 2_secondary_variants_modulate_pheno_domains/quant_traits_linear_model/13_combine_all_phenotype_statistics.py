#!/bin/python3



import pandas as pd


df = pd.DataFrame()


# Pearson's R
filename = 'statistics/correlation_phenotypes.csv'
app = pd.read_csv(filename)
app.columns = ['Variant class/PRS', 'Phenotype', 'Pvalue', 'Sample size', 'Direction', 'FDR']
app['Direction metric'] = 'Pearson\'s R'
app['Test'] = 'Pearson\'s R'
app['Other info'] = ''
app['Group'] = '16p12_WGS'
df = df.append(app)

# T-test
filename = 'statistics/ttest_phenotypes.csv'
app = pd.read_csv(filename)
app.columns = ['Variant class/PRS', 'Phenotype', 'Pvalue', 'Sample size', 'Direction', 'FDR']
app['Direction metric'] = 'Cohen\'s D'
app['Test'] = 'T-test'
app['Other info'] = ''
app['Group'] = '16p12_WGS'
df = df.append(app)

# T-test Family history
filename = 'statistics/ttest_family_history.csv'
app = pd.read_csv(filename)
app.columns = ['Variant class/PRS', 'Phenotype', 'Pvalue', 'Sample size', 'Direction', 'FDR']
app['Direction metric'] = 'Cohen\'s D'
app['Test'] = 'T-test'
app['Other info'] = ''
app['Group'] = '16p12_WGS'
df = df.append(app)

# Anova
filename = 'statistics/anova_phenotypes.csv'
app = pd.read_csv(filename)
app.columns = ['Variant class/PRS', 'Phenotype', 'Pvalue', 'Sample size', 'Direction', 'FDR']
app['Direction metric'] = 'Partial omega squared'
app['Test'] = 'ANNOVA'
app['Other info'] = ''
app['Group'] = '16p12_WGS'
df = df.append(app)



       
outfile = 'statistics/all_statistics_phenotypes.xlsx'
df.to_excel(outfile, index=False)


