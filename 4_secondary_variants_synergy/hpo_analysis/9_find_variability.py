import pandas as pd

# Find the number of genes in our cohort that are associated with different phenotypes in differnet probands

# Load in files
# Phenotypes
pheno = pd.read_csv('Annotated_files/phenotype_hpo.csv')
# HPO ID - phenotype map
mc2hpo = pd.read_excel('Analysis_files/MC_HPO_terms.xlsx')
mc2hpo.dropna(how = 'any', inplace = True)
mc2hpo = mc2hpo[['HPO ID', 'HPO Term']]
mc2hpo.drop_duplicates(inplace = True, keep = 'first')
mc2hpo.set_index(mc2hpo['HPO ID'], inplace = True)
hpo_dict = mc2hpo['HPO Term'].to_dict()
# Variants
snvs = pd.read_csv('Annotated_files/snvs_hpo_anno.csv')
cnvs = pd.read_csv('Annotated_files/cnvs_hpo_anno.csv')
strs = pd.read_csv('Annotated_files/strs_hpo_anno.csv')
# Make gene symbol, gene_id map
df = cnvs[['Sample', 'Gene_Symbol', 'Gene_ID', 'GeneHPO_inclusive']]
df.columns = ['Sample', 'Gene_symbol', 'Gene_id_', 'GeneHPO_inclusive']
df = pd.concat([snvs[['Sample', 'Gene_symbol', 'Gene_id_', 'GeneHPO_inclusive']], df, strs[['Sample', 'Gene_symbol', 'Gene_id_', 'GeneHPO_inclusive']]])

# Filter for only probands
cohort_info=pd.read_csv('../../11_Variant Integration/16p12_cohort_summary_v21.csv')
cohort_info=cohort_info[cohort_info.No_consent_forms.isnull()]
cohort_info.index=cohort_info.Sample
df['Relationship']=df.Sample.map(cohort_info.Relationship.to_dict())

df=df[df.Relationship=='P']

# Filter for probands with phenotypes
df=df[df.Sample.isin(pheno.Sample.to_list())]

# Restrict to recurrent genes
recurr_genes=df.Gene_symbol.value_counts()
recurr_genes=recurr_genes[recurr_genes>=2]

df=df[df.Gene_symbol.isin(recurr_genes.index.to_list())]
df.sort_values('Gene_symbol', inplace=True)

# Restrict to genes with HPO annotations
df=df[df.GeneHPO_inclusive!='.']

# Get associated phenotypes from proband
def get_associated(row):
    gene_terms = row.GeneHPO_inclusive
    gene_hpo = row.GeneHPO_inclusive.split(';')

    pro_terms = pheno[pheno.Sample==row.Sample]['HPO_ID'].to_list()[0].split(';')

    overlap_terms = list(set([i for i in pro_terms if i in gene_hpo]))
    if overlap_terms==[]:
        overlap_terms=['.']
    overlap_terms.sort()
    
    return ';'.join(overlap_terms)
df['Associated_phenotypes']=df.apply(get_associated, axis=1)

# Get number of genes with different associated phenotypes in each proband
phenotype_genes=list(set(df.Gene_symbol.to_list()))
print('Total number of recurrent genes:', len(phenotype_genes))
count=0
for g in phenotype_genes:
    pheno_assoc=list(set(df[df.Gene_symbol==g]['Associated_phenotypes'].to_list()))
    if len(pheno_assoc)>1:
        count+=1
print('Total number of recurrent genes with variable phenotypes:', str(count))

# Restrict to probands with at least 3 phenotypes
pheno3=pheno[pheno.number_of_phenotypes>=3]
df=df[df.Sample.isin(pheno3.Sample.to_list())]

# Restrict to recurrent genes
recurr_genes=df.Gene_symbol.value_counts()
recurr_genes=recurr_genes[recurr_genes>=2]

df=df[df.Gene_symbol.isin(recurr_genes.index.to_list())]
df.sort_values('Gene_symbol', inplace=True)

phenotype_genes=list(set(df.Gene_symbol.to_list()))
print('Total number of recurrent genes in probands with >=3 phenotypes:', len(phenotype_genes))
count=0
for g in phenotype_genes:
    pheno_assoc=list(set(df[df.Gene_symbol==g]['Associated_phenotypes'].to_list()))
    if len(pheno_assoc)>1:
        count+=1
print('Total number of recurrent genes with variable phenotypes in probands with >=3 phenotypes:', str(count))