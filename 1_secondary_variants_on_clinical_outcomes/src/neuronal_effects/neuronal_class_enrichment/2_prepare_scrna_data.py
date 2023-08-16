#!/bin/python3


import pandas as pd




cell_df = pd.read_csv('allen_institute_tables/metadata.csv')
df = pd.read_csv('allen_institute_tables/medians.csv')


# the column cluster_label and subclass_label
# seem to hold the information about the cell type
cell_df = cell_df[['cluster_label', 'subclass_label']]
cell_df['superclass_label'] = cell_df["cluster_label"].str.split(" ", 1).str[0]
cell_df = cell_df.drop_duplicates()
cell_df = cell_df.set_index('cluster_label', drop=True)

#Save file that maps the subclasses to the superclasses for future reference
cell_class_df=cell_df
cell_class_df = cell_class_df.set_index('superclass_label',drop=True)
cell_class_df = cell_class_df.drop_duplicates()
cell_class_df.to_csv('allen_institute_tables/cell_class_mapping.csv', index=True)

#Superclass label = 8 classes; Subclass label = 20 classes


# subclass_label
cell_type_labels = list(df.columns[1:])
# add subclass information to the expression matrix column names
#cell_type_labels = [s+'|'+cell_df.at[s, 'subclass_label'] for s in cell_type_labels]
df = df.set_index('feature', drop=True)
df.columns = cell_type_labels
#Add subclass and superclass labels
df=pd.concat([df,cell_df[["subclass_label"]].transpose()])
df=pd.concat([df,cell_df[["superclass_label"]].transpose()])

superclass_list=list(set(cell_df['superclass_label'].values.tolist()))
subclass_list=list(set(cell_df['subclass_label'].values.tolist()))

#Generate summed read counts for each subclass type

subclass_df=pd.DataFrame()
for subclass in subclass_list:
	sub_df = df.loc[:, df.loc['subclass_label'] == subclass]
	sub_df = sub_df.drop(['superclass_label','subclass_label'])
	gene_sum = sub_df.sum(axis=1)
	subclass_df[subclass]=gene_sum

# median expression of each gene
# and standard deviation for each gene
gene_medians = subclass_df.median(axis=1)
gene_std_dev = subclass_df.std(axis=1)


# create a new df that will hold cell type specific expression information
# A value of 1 means that gene is preferentially expressed in the cell type
# a value of 0 means that it is not preferentially expressed in the cell type
new_df = pd.DataFrame(index=subclass_df.index)
for cell_type in subclass_df.columns:
	new_df[cell_type] = subclass_df[cell_type] > (gene_medians + 2* gene_std_dev)
	new_df[cell_type] = new_df[cell_type].astype(int)
	print(new_df[cell_type].sum(), cell_type) # print the number of pref exp. genes for this cell type

new_df.to_csv('allen_institute_tables/cell_type_specific_expression_subclass.csv', index=True)



#Repeat above for superclasses of cell types

superclass_df=pd.DataFrame()
for superclass in superclass_list:
	subdf = df.loc[:, df.loc['superclass_label'] == superclass]
	subdf = subdf.drop(['superclass_label','subclass_label'])
	gene_sum = subdf.sum(axis=1)
	superclass_df[superclass]=gene_sum

# median expression of each gene
# and standard deviation for each gene
gene_medians = superclass_df.median(axis=1)
gene_std_dev = superclass_df.std(axis=1)


# create a new df that will hold cell type specific expression information
# A value of 1 means that gene is preferentially expressed in the cell type
# a value of 0 means that it is not preferentially expressed in the cell type
new_df = pd.DataFrame(index=superclass_df.index)
for cell_type in superclass_df.columns:
	new_df[cell_type] = superclass_df[cell_type] > (gene_medians + 2* gene_std_dev)
	new_df[cell_type] = new_df[cell_type].astype(int)
	print(new_df[cell_type].sum(), cell_type) # print the number of pref exp. genes for this cell type

new_df.to_csv('allen_institute_tables/cell_type_specific_expression_superclass.csv', index=True)



