import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['pdf.fonttype'] = 42

df=pd.read_csv('result_tables/1_covariate_even_split.csv')
# Remove variables we don't want to plot
df=df[df.Test!='Linear regression']
df=df[~df.Variable.isin(['(Intercept)', 'Sex'])]
df=df[df.model!='Sex + All_rare_del_var + SCZ_PRS']
df=df[(df.model.str.contains('SCZ_PRS')) & (df.model!='Sex + SCZ_PRS')]

df['num_vars']=df.model.apply(lambda m: len(m.split('+')))

df['model_short']='Model'
df.loc[df.model.str.contains('LOEUF'), 'model_short']='Model (LF)'
df.loc[df.num_vars==2, 'model_short']='Single'

df=df[df.model_short!='Single']

domains={'Child_behav':'Behavioral', 'Child_psych':'Psychiatric', 'Child_nervous_system':'Nervous system', 'Child_congenital':'Congenital', 'Child_craniofacial':'Growth/Skeletal'}
df['domain']=pd.Categorical(df.phenotype.map(domains), ['Behavioral', 'Psychiatric', 'Nervous system', 'Congenital', 'Growth/Skeletal'])
domain_colors={'Behavioral':'#9586BF', 'Psychiatric':'#FAB263', 'Nervous system':'#82B2D3', 'Congenital':'#FAF6B7', 'Growth/Skeletal':'#F17F72'}
vars={'Rare_Deleterious_SNVs':'All coding SNVs', 'genes_del':'Genes del.', 'genes_dup':"Genes dup", 'STRs_exonic':'STRs',
		'Rare_Deleterious_SNVs_LOEUF':'All coding SNVs (LF)', 'dels_loeuf':'Genes del. (LF)', 'dups_loeuf':'Genes dup. (LF)', 'STRs_exonic_LOEUF035':'STRs (LF)', 'SCZ_PRS':'SCZ PRS'}
df['var']=pd.Categorical(df.Variable.map(vars), ['All coding SNVs', 'Genes del.', "Genes dup", 'STRs', 'All coding SNVs (LF)', 'Genes del. (LF)', 'Genes dup. (LF)', 'STRs (LF)', 'SCZ PRS'])

# Add some dummay variables for later
for v in ['All coding SNVs', 'Genes del.', "Genes dup", 'STRs', 'SCZ PRS']:
	df=pd.concat([df, pd.DataFrame({'var':[v]*2, 'model_short':['Model']*2, 'domain':['Z']*2})])
for v in ['All coding SNVs (LF)', 'Genes del. (LF)', 'Genes dup. (LF)', 'STRs (LF)', 'SCZ PRS']:
	df=pd.concat([df, pd.DataFrame({'var':[v]*2, 'model_short':['Model (LF)']*2, 'domain':['Z']*2})])

df.sort_values(by=['model_short', 'var', 'domain'], inplace=True, ascending=[True, False, True])

pdf=PdfPages('Figures/3_forest_plots.pdf')
for mod in ['Model', 'Model (LF)']:
	fig, ax=plt.subplots(figsize=(6, 8))
	
	subdf=df[df.model_short==mod].copy()
	subdf.reset_index(drop=True, inplace=True)
	subdf['y']=subdf.index.to_list()
	subdf=subdf[subdf.domain!='Z']
	
	print(subdf)
	
	for d in ['Behavioral', 'Psychiatric', 'Nervous system', 'Congenital', 'Growth/Skeletal']:
		subdf2=subdf[subdf.domain==d].copy()
		# Ranges
		for i in range(5):
			plt.plot((subdf2['2.5% C.I.'].to_list()[i], subdf2['97.5% C.I.'].to_list()[i]), (subdf2['y'].to_list()[i], subdf2['y'].to_list()[i]), color=domain_colors[d])
		# Center dots
		plt.scatter(subdf2['regression_coefficient'].to_list(), subdf2['y'].to_list(), color=domain_colors[d])
	
	# Add a line at 0
	lo, hi=plt.ylim()
	plt.plot((0, 0), (lo-1, hi+1), color='k', ls=':')
	plt.ylim(lo, hi)
	
	# Add axis labels
	if mod=='Model':
		plt.yticks([30, 23, 16, 9, 2], ['All coding SNVs', 'Genes del.', "Genes dup", 'STRs', 'SCZ PRS'])
	else:
		plt.yticks([30, 23, 16, 9, 2], ['All coding SNVs (LF)', 'Genes del. (LF)', 'Genes dup. (LF)', 'STRs (LF)', 'SCZ PRS'])
	
	# Legend
	custom_lines=[]
	for d in ['Behavioral', 'Psychiatric', 'Nervous system', 'Congenital', 'Growth/Skeletal']:
		custom_lines.append(Line2D([0], [0], color=domain_colors[d], lw=4))
	plt.legend(custom_lines, ['Behavioral', 'Psychiatric', 'Nervous system', 'Congenital', 'Growth/Skeletal'], bbox_to_anchor=(1, 1))
	
	# Add title and format layout
	plt.title(mod)
	plt.tight_layout()

	pdf.savefig()
	plt.close()
pdf.close()
		