import pandas as pd
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

matplotlib.rcParams['pdf.fonttype'] = 42

pdf=PdfPages('Figures/2_model_results.pdf')

for split in ['uneven', 'even']:
    df=pd.read_csv('result_tables/1_covariate_'+split+'_split.csv')
    
    # Remove variables we don't want to plot
    df=df[~df.Variable.isin(['(Intercept)', 'Sex'])]
    df=df[df.model!='Sex + All_rare_del_var + SCZ_PRS']
    df=df[(df.model.str.contains('SCZ_PRS')) & (df.model!='Sex + SCZ_PRS')]
    df=df[~((df.Variable=='SCZ_PRS') & (df.model=='Sex + SCZ_PRS + Rare_Deleterious_SNVs_LOEUF + dels_loeuf + dups_loeuf + STRs_exonic_LOEUF035'))]
    
    for mod in ['Logistic regression', 'Linear regression']:
        subdf=df[df.Test==mod]
        
        vars=['Rare_Deleterious_SNVs', 'genes_del', 'genes_dup', 'STRs_exonic', 'SCZ_PRS', 'Rare_Deleterious_SNVs_LOEUF', 'dels_loeuf', 'dups_loeuf', 'STRs_exonic_LOEUF035']
        phenos=['Child_behav', 'Child_psych', 'Child_nervous_system', 'Child_congenital', 'Child_craniofacial']
        if mod=='Linear regression':
            phenos=['De_vrie', 'srs', 'hrs_mat', 'head_circumference_z_score', 'bmi_z_score']
        
        # Add p-value annotations
        subdf['p_anno']=''
        subdf.loc[subdf['P-value']<=0.05, 'p_anno']='*'
        
        # Convert to wide-form data
        wide_vals=subdf.pivot(index='Variable', columns='phenotype', values='regression_coefficient')
        wide_vals=wide_vals.loc[vars, phenos]
        wide_ps=subdf.pivot(index='Variable', columns='phenotype', values='p_anno')
        wide_ps=wide_ps.loc[vars, phenos]
        
        # Draw heatmap
        colors=["#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"]
        cmap=LinearSegmentedColormap.from_list('Blue-Red', colors, N=20)
        sns.heatmap(data=wide_vals, cmap=cmap, vmin=-1, vmax=1, center=0, annot=wide_ps, square=True, fmt='', linecolor='grey', linewidths=0.75)
        
        plt.title(split+'_'+mod)
        plt.tight_layout()
        pdf.savefig()
        plt.close()

pdf.close()

# Make forest plots for each equation and phenotype
def forest(subdf, title):
    fig, ax=plt.subplots(figsize=(4, 5))
    ys=[i for i in range(subdf.shape[0])]
    ys.reverse()
    est=subdf.regression_coefficient.to_list()
    lower=subdf['2.5% C.I.'].to_list()
    upper=subdf['97.5% C.I.'].to_list()
    
    plt.scatter(est, ys, color='k')
    for i in range(len(ys)):
        plt.plot((lower[i], upper[i]), (ys[i], ys[i]), color='k')
            
    # Add a line at 0
    lo, hi=plt.ylim()
    plt.plot((0, 0), (lo-1, hi+1), color='k', ls=':')
    plt.ylim(lo, hi)
    
    # Add axis labels
    plt.yticks(ys, subdf.Variable.to_list())
    
    # Add title and format layout
    plt.title(title)
    plt.tight_layout()

for split in ['uneven', 'even']:
    df=pd.read_csv('result_tables/1_covariate_'+split+'_split.csv')
    
    # Remove variables we don't want to plot
    df=df[df.Variable!='(Intercept)']
    df=df[df.model!='Sex + All_rare_del_var + SCZ_PRS']
    df=df[(df.model.str.contains('SCZ_PRS')) & (df.model!='Sex + SCZ_PRS')]
    
    # Make variables categorical for easy sorting
    df.Variable=pd.Categorical(df.Variable, ['Sex', 'Rare_Deleterious_SNVs', 'genes_del', 'genes_dup', 'STRs_exonic', 'Rare_Deleterious_SNVs_LOEUF', 'dels_loeuf', 'dups_loeuf', 'STRs_exonic_LOEUF035', 'SCZ_PRS'])
    
    pdf=PdfPages('Figures/2_model_forest_'+split+'.pdf')
    
    # Make plots with and without covariates
    for mod in list(df.model.unique()):
        name='model'
        if 'LOEUF' in mod:
            name='model (LF)'
        
        for pheno_lst in [['Child_behav', 'Child_psych', 'Child_nervous_system', 'Child_congenital', 'Child_craniofacial'], ['De_vrie', 'srs', 'hrs_mat', 'head_circumference_z_score', 'bmi_z_score']]:
            for p in pheno_lst:
                subdf=df[(df.model==mod) & (df.phenotype==p)]
                subdf.sort_values('Variable', inplace=True)
            
                forest(subdf, name+' - '+p)
                pdf.savefig()
                plt.close()
    
                subdf=subdf[~subdf.Variable.isin(['Sex'])]
                forest(subdf, name+' - '+p+' - no sex')
                pdf.savefig()
                plt.close()

    pdf.close()
    
# Plot variance
pdf=PdfPages('Figures/2_model_variance.pdf')
for split in ['even', 'uneven']:
    df=pd.read_csv('result_tables/1_covariate_'+split+'_split.csv')
    
    # Remove variables we don't want to plot
    df=df[~df.Variable.isin(['(Intercept)', 'Sex'])]
    df=df[df.model!='Sex + All_rare_del_var + SCZ_PRS']
    df=df[~((df.Variable=='SCZ_PRS') & (df.model=='Sex + SCZ_PRS + Rare_Deleterious_SNVs_LOEUF + dels_loeuf + dups_loeuf + STRs_exonic_LOEUF035'))]
    
    # Add in shortened names
    name_dict={'Sex + All_rare_del_var':'Model1',
                'Sex + SCZ_PRS + Rare_Deleterious_SNVs + genes_del + genes_dup + STRs_exonic':'All coding variants',
                'Sex + Rare_Deleterious_SNVs':'All coding SNVs',
                'Sex + genes_del':'Genes del.',
                'Sex + genes_dup':'Genes dup',
                'Sex + STRs_exonic':'STRs',
                'Sex + SCZ_PRS':'SCZ PRS',
                'Sex + SCZ_PRS + Rare_Deleterious_SNVs_LOEUF + dels_loeuf + dups_loeuf + STRs_exonic_LOEUF035':'All coding variants (LF)',
                'Sex + Rare_Deleterious_SNVs_LOEUF':'All coding SNVs (LF)',
                'Sex + dels_loeuf':'Genes del. (LF)',
                'Sex + dups_loeuf':'Genes dup. (LF)',
                'Sex + STRs_exonic_LOEUF035':'STRs (LF)'}
    order=[name_dict[i] for i in name_dict.keys()]
    df['name']=df.model.map(name_dict)
    
    for mod in ['Logistic regression', 'Linear regression']:
        subdf=df[df.Test==mod]
        
        phenos=['Child_behav', 'Child_psych', 'Child_nervous_system', 'Child_congenital', 'Child_craniofacial']
        if mod=='Linear regression':
            phenos=['De_vrie', 'srs', 'hrs_mat', 'head_circumference_z_score', 'bmi_z_score']
        
        # Limit to relevant entries
        subdf=subdf[['name', 'model', 'phenotype', 'R2']]
        subdf.drop_duplicates(inplace=True)
        
        # Convert to wide-form data
        wide_vals=subdf.pivot(index='name', columns='phenotype', values='R2')
        wide_vals=wide_vals.loc[order, phenos]
        
        # Draw heatmap
        colors=["#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"]
        cmap=LinearSegmentedColormap.from_list('Reds', colors, N=20)
        sns.heatmap(data=wide_vals, cmap=cmap, vmin=0, vmax=0.2, square=True, fmt='', linecolor='grey', linewidths=0.75)
        
        plt.title(split+'_'+mod)
        plt.tight_layout()
        pdf.savefig()
        plt.close()
pdf.close()