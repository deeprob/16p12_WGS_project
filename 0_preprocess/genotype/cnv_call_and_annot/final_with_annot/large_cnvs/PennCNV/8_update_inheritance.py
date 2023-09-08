import pandas as pd

# Using the hand annotations from the visualizations (see script 7), update the inheritance annotations for the calls
calls = pd.read_csv('call_tables/6_genic_filter.txt', sep = '\t')
put_dn = pd.read_csv('call_tables/7.1_denovo_anno.txt', sep = '\t')

df = pd.merge(calls, put_dn, on = calls.columns.to_list(), how = 'left')

print(df)

# Update inheritance
def update_inherit(inherit1, inherit2):
	if inherit1=='both':
		return '.'
	elif inherit1!='de novo':
		return inherit1
	else:
		return inherit2

df['inheritance'] = df[['inheritance', 'Actual_inheritance']].apply(lambda row: update_inherit(row[0], row[1]), axis = 1)
print(df)

df.to_csv('call_tables/8_update_inheritance.txt', sep = '\t', index = False)
