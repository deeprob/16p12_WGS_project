import pandas as pd

# Make a table with the following information:
# Sample, Gene, Mutation type, HPO terms associated with gene, Percentage phenotypes explained by gene

# Get patient phenotypes
pheno=pd.read_csv('Annotated_files/phenotype_hpo.csv')

# Get variants and associated HPO terms
def number_anno(row):
    samp=row['Sample']
    gene_terms=row['GeneHPO_inclusive'].split(';')
    
    if samp not in pheno.Sample.to_list():
        return 'No available phenotypes in patient'
    samp_terms=pheno[pheno.Sample==samp]['HPO_ID'].to_list()[0].split(';')
    if samp_terms==['0']:
        return 'No phenotypes in patient'
    explained_terms=[i for i in samp_terms if i in gene_terms]
    return len(explained_terms)/len(samp_terms)

vars=['snvs', 'cnvs', 'strs']
outdf=pd.DataFrame(0, index=[0], columns=['Sample', 'Gene_id_', 'Gene_symbol', 'Mut_type', 'GeneHPO_inclusive', 'Percentage_explained'])
for v in vars:
    df=pd.read_csv('Annotated_files/'+v+'_hpo_anno.csv')
    if v=='cnvs':
        df['Mut_type']=df.Type
        df['Gene_id_']=df.Gene_ID
        df['Gene_symbol']=df.Gene_Symbol
    elif v=='strs':
        df['Mut_type']='STR'
    df=df[['Sample', 'Gene_id_', 'Gene_symbol', 'Mut_type', 'GeneHPO_inclusive']]
    
    df['Percentage_explained']=df.apply(number_anno, axis=1)
    outdf=pd.concat([outdf, df])

# Clean up and save to file
outdf=outdf[~outdf.Percentage_explained.isin(['No available phenotypes in patient', 'No phenotypes in patient'])]
outdf=outdf[outdf.Sample!=0]
outdf.loc[outdf.GeneHPO_inclusive=='.', 'GeneHPO_inclusive']='No HPO terms associated with gene'
outdf.sort_values(by=['Sample', 'Mut_type', 'Gene_id_'], inplace=True)
outdf.columns=['Sample', 'Ensembl ID', 'Gene symbol', 'Mutation type', 'HPO terms associated with gene', 'Percentage phenotypes explained by gene']
outdf.to_csv('Result_tables/8_tableS3.csv', index=False)