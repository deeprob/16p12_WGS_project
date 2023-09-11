import pandas as pd
import numpy as np

# Group probands 2 ways:
# 1. No parent phenotypes, One parent has phenotypes, 2 parents have phenotypes

df=pd.read_csv('Analysis_files/1_domain_scores.csv')

m_cols=[i for i in df.columns.to_list() if 'motherParent'in i]
f_cols=[i for i in df.columns.to_list() if 'fatherParent'in i]
keep_cols=['Proband', 'Family', 'Mother', 'Father']+m_cols+f_cols
df=df[keep_cols]

# Drop families missing any values
df.dropna(how='any', axis=0, inplace=True)

df['mother_sum']=df[m_cols].sum(axis=1)
df['father_sum']=df[f_cols].sum(axis=1)

df.loc[(df.mother_sum==0) & (df.father_sum==0), 'number_parents']='C'
df.loc[(df.mother_sum>0) | (df.father_sum>0), 'number_parents']='B'
df.loc[(df.mother_sum>0) & (df.father_sum>0), 'number_parents']='A'

print(df.number_parents.value_counts())

# 2. No phenotypes, each phenotypic domain, and comorbid domains
domains=['SCZ', 'depression', 'cognitive', 'addiction']
group_cols=[]
for d in domains:
    df['parent_sum_'+d] = df[['motherParent_'+d, 'fatherParent_'+d]].sum(axis=1)
    df.loc[df['parent_sum_'+d]==0, 'parent_has_'+d] = False
    df.loc[df['parent_sum_'+d]>0, 'parent_has_'+d] = True
    df.drop('parent_sum_'+d, inplace=True, axis=1)
    group_cols.append('parent_has_'+d)
count_df=df.groupby(group_cols).size()
print(count_df)

# Based on sizes of groups, breakdown will be: No parent phenotypes, No SCZ phenotype, Other phenotypes
df.loc[~(df.parent_has_cognitive & df.parent_has_depression & df.parent_has_SCZ & df.parent_has_addiction), 'phenotype_parents']='C'
df.loc[~(df.parent_has_SCZ) & (df.parent_has_depression | df.parent_has_cognitive | df.parent_has_addiction), 'phenotype_parents']='B'
df.loc[df.parent_has_SCZ, 'phenotype_parents']='A'
print(df.phenotype_parents.value_counts())

# Save clusters to file
df[['Proband']+m_cols+f_cols+['number_parents', 'phenotype_parents']].to_csv('Analysis_files/3_grouped_clusters.csv', index=False)