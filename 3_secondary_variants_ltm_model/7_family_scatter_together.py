import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
sns.set_style({'font.family':'sans-serif', 'font.sans-serif':'Arial'})
from statannot import add_stat_annotation
from matplotlib.backends.backend_pdf import PdfPages

# Get correlation statistics
stats=pd.read_csv('../../20_Burden-correlation analysis/Correlation_Matrices/probands_and_parents/statistics/variant_v_variant.csv')
stats=stats[~stats.FDR.isnull()]
stats=stats[(stats.Variant1.str.contains('PRS')) & (stats.Variant2=='Rare_Deleterious_SNVs')]

# Get raw data
df=pd.read_csv('../../11_Variant Integration/16p12_cohort_summary_v22.csv')
df['rel']='.'
df.loc[df.Relationship=='P', 'rel']='Proband'
df.loc[df.Relationship.isin(['MC', 'FC']), 'rel']='Carrier Parent'
df.loc[df.Relationship.isin(['MNC', 'FNC']), 'rel']='Noncarrier Parent'
df=df[df.rel!='.']
df=df[~df.Rare_Deleterious_SNVs.isnull()]

pdf=PdfPages('figures/7_correlation_together.pdf')

prs=list(stats.Variant1.unique())
groups=['Carrier Parent', 'Noncarrier Parent']
group_colors={'Carrier Parent':'#991323', 'Noncarrier Parent':'#26a9e0'}

for p in prs:
    fig, axs = plt.subplots(nrows=1,ncols=2, figsize=(8,3), sharex=True, sharey=True)
    
    for i in range(2):
        g=groups[i]
        subdf=df[(~df[p].isnull()) & (df.rel.isin([g, 'Proband']))]
        
        # Plot scatter and best fit line
        sns.scatterplot(data=subdf, x=p, y='Rare_Deleterious_SNVs', hue='rel', hue_order=['Proband', g], palette=['#2a8057', group_colors[g]], ax=axs[i])
        
        sns.regplot(data=subdf[subdf.rel=='Proband'], x=p, y='Rare_Deleterious_SNVs', color='#2a8057', ax=axs[i], scatter=False)
        sns.regplot(data=subdf[subdf.rel==g], x=p, y='Rare_Deleterious_SNVs', color=group_colors[g], ax=axs[i], scatter=False)
        
        # Annotate stats
        for fg in ['Proband', g]:
            substat=stats[(stats.Variant1==p) & (stats.Sample_Group==fg)]
            r=substat.Coeff.to_list()[0]
            pval=substat.Pvalue.to_list()[0]
            fdr=substat.FDR.to_list()[0]
            text=('R=%.2f\np=%.1E\nFDR=%.2f' % (r, pval, fdr))
            if fg=='Proband':
                axs[i].annotate(text, xy=(0.1, 0.75), xycoords='axes fraction')
            else:
                axs[i].annotate(text, xy=(0.7, 0.75), xycoords='axes fraction')
        axs[i].legend(bbox_to_anchor=(1, 1))
        
        axs[i].set_title(g)
    pdf.savefig()
    plt.close()

pdf.close()