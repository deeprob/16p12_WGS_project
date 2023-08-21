import pandas as pd

from upsetplot import plot
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42

# Create an upset plot showing the probands with mutations in each category of Pathogenic variant
# Categories are:
# 1. LOF mutations in Tier S or 1 SFARI genes
# 2. LOF mutations in Tier 1 or 2 DBD genes
# 3. LOF mutations in genes with dominant or X-linked OMIM annotations with relevant phenotypes (from a manual screen)
# 4. "Pathogenic" or "Likely pathogenic" variant in ClinVar with relevant phenotypes/mode of action (from a manual screen)

df=pd.read_csv('Rare_Deleterious_Exonic_Variants_Pathogenic_Anno.csv')

df.replace('.', 0, inplace=True)
df.replace('X', 1, inplace=True)

# Restrict to probands
cohort_info=pd.read_csv('../../11_Variant Integration/16p12_cohort_summary_v21.csv')
cohort_info=cohort_info[cohort_info.No_consent_forms.isnull()]
cohort_info.index=cohort_info.Sample
df['Relationship']=df.Sample.map(cohort_info.Relationship.to_dict())

df=df[df.Relationship=='P']

# Add in NEJM CNV data
cnv=pd.read_csv('../../9_Rare variant calls/Rare coding CNV calls/sv_calls_combined.txt', sep='\t')
cnv=cnv[cnv.NEJM!='.'][['Sample', 'NEJM']]
cnv=cnv[cnv.Sample.isin(df.Sample.to_list())]
cnv=cnv[~((cnv.NEJM.str.contains('p12.1')) & (cnv.NEJM.str.contains('16p')))]
cnv.drop_duplicates(inplace=True)

df['Pathogenic_CNV']=0
df.loc[df.Sample.isin(cnv.Sample.to_list()), 'Pathogenic_CNV']=1

# Condense data on the proband level
proband_df=df[['Sample', 'ClinVar_Manual_Screen', 'OMIM_Manual_Screen', 'SFARI_Tier_1_S', 'DBD_Tier_1_2', 'Pathogenic_CNV']].groupby('Sample').sum()
for col in ['ClinVar_Manual_Screen', 'OMIM_Manual_Screen', 'SFARI_Tier_1_S', 'DBD_Tier_1_2', 'Pathogenic_CNV']:
    proband_df.loc[proband_df[col]>0, col]=True
    proband_df.loc[proband_df[col]==0, col]=False
proband_df['None']=True
proband_df.loc[(proband_df.ClinVar_Manual_Screen) | (proband_df.OMIM_Manual_Screen) | (proband_df.SFARI_Tier_1_S) | (proband_df.DBD_Tier_1_2) | (proband_df.Pathogenic_CNV), 'None']=False

# Update column names for plot
proband_df=proband_df[['ClinVar_Manual_Screen', 'OMIM_Manual_Screen', 'SFARI_Tier_1_S', 'DBD_Tier_1_2', 'Pathogenic_CNV', 'None']]
proband_df.columns=['ClinVar', 'OMIM (Dominant/X-linked)', 'Tier S/1 SFARI', 'Tier 1/2 DBD', 'Pathogenic CNV', 'None']

# Make upset plot
plot_df=proband_df.groupby(['ClinVar', 'OMIM (Dominant/X-linked)', 'Tier S/1 SFARI', 'Tier 1/2 DBD', 'Pathogenic CNV', 'None']).size()
print(plot_df)
plot(plot_df)
plt.savefig('Results/3_pathogenic_variant_upset.pdf')

# Save data to file
plot_df.to_csv('Results/3_upset_data.csv')