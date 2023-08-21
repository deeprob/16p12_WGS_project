import pandas as pd

# Get the phenotypes associated with the ClinVar variant

# Load in data
# Get probands
cohort_info = pd.read_csv('../../3_cohort information/16p12_All_Participants_v9.csv')
cohort_info = cohort_info[(cohort_info.WGS=='X') & (cohort_info.No_consent_forms!='X')]
pros = cohort_info[cohort_info.Relationship=='P']['Sample'].to_list()

# Phenotypes
pheno = pd.read_csv('../1_HPO_annotations/Annotated_files/phenotype_hpo.csv')[['Sample', 'HPO_ID']]
pheno = pheno[pheno.Sample.isin(pros)]
# Variants
snvs = pd.read_csv('Analysis_files/4_ClinVar_variants_manualHPO.csv')[['Sample', 'Gene_symbol', 'Mut_type', 'AAChange.wgEncodeGencodeBasicV19', 'ClinVar_ALLELEID', 'ClinVar_Manual_Screen', 'Gene_id_', 'VarHPO', 'GeneHPO_inclusive']]
snvs = snvs[(snvs.Sample.isin(pros)) & (snvs.ClinVar_Manual_Screen=='X')]
snvs.sort_values('Sample', inplace=True)
print(snvs)

pheno=pheno[pheno.Sample.isin(snvs.Sample.to_list())]

# Get HPO translations
mc2hpo = pd.read_excel('../1_HPO_annotations/Analysis_files/MC_HPO_terms.xlsx')
mc2hpo.dropna(how = 'any', inplace = True)
mc2hpo = mc2hpo[['HPO ID', 'HPO Term', 'Phenotype Category']]
mc2hpo.drop_duplicates(inplace = True, keep = 'first')
mc2hpo.set_index(mc2hpo['HPO ID'], inplace = True)
hpo_dict = mc2hpo['HPO Term'].to_dict()

# Get groups of related phenotypes
domain_defs = dict(zip(mc2hpo['HPO ID'].to_list(), mc2hpo['Phenotype Category'].to_list()))
domain_defs2=dict(zip(mc2hpo['HPO Term'].to_list(), mc2hpo['Phenotype Category'].to_list()))
domains = list(set(mc2hpo['Phenotype Category'].to_list()))

def get_pheno(pro, gene):
    pro_pheno=pheno[pheno.Sample==pro]['HPO_ID'].to_list()[0].split(';')
    pro_pheno_sorted=[]
    for d in domains:
        domain_phenos=[i for i in pro_pheno if domain_defs[i]==d]
        domain_phenos.sort()
        pro_pheno_sorted+=domain_phenos
    
    var_pheno=snvs[(snvs.Sample==pro) & (snvs['Gene_symbol']==gene)]['VarHPO'].to_list()[0].split(';')
    
    overlap_ids=list(set(pro_pheno_sorted).intersection(set(var_pheno)))
    overlap_phenos=[hpo_dict[i] for i in overlap_ids]
    
    proband_only_ids=list(set(pro_pheno_sorted)-set(var_pheno))
    proband_only_phenos=[hpo_dict[i] for i in proband_only_ids]
    
    print('Proband phenotypes associated with variant:', ', '.join(overlap_phenos))
    print('Proband phenotypes not associated with variant:', ', '.join(proband_only_phenos))
    
    overlap_domains=list(set([domain_defs2[i] for i in overlap_phenos]))
    overlap_domains.sort()
    
    proband_only_domains=list(set([domain_defs2[i] for i in proband_only_phenos]))
    proband_only_domains.sort()
    
    # Also print domains
    print('Proband domains associated with variant:', ', '.join(overlap_domains))
    print('Proband phenotypes not associated with variant:', ', '.join(proband_only_domains))

# Get phenotypes associated with variants in each sample
samps=list(snvs.Sample.unique())
for samp in samps:
    genes=snvs[snvs.Sample==samp]['Gene_symbol'].to_list()
    for g in genes:
        cvai=str(snvs[(snvs.Sample==samp) & (snvs.Gene_symbol==g)]['ClinVar_ALLELEID'].to_list()[0])+'[AlleleID]'
        print(samp, g, cvai) 
        get_pheno(samp, g)
        print('\n')
    
        # Print domains in each child
        
    