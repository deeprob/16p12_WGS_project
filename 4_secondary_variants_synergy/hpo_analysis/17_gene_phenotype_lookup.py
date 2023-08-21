#!/bin/python
import pandas as pd

from thefuzz import fuzz

# A script to quickly lookup the phenotypes explained by a gene for a proband
# This script will be used to manually annotate the figures in Adobe

# Load in files
# Phenotypes
pheno = pd.read_csv('../1_HPO_annotations/Annotated_files/phenotype_hpo.csv')[['Sample', 'HPO_ID']]
# HPO ID - phenotype map
mc2hpo = pd.read_excel('../1_HPO_annotations/Analysis_files/MC_HPO_terms.xlsx')
mc2hpo.dropna(how = 'any', inplace = True)
mc2hpo = mc2hpo[['MC Phenotype', 'HPO ID', 'HPO Term']]
mc2hpo.drop_duplicates(inplace = True, keep = 'first')
mc2hpo.index = mc2hpo['HPO ID'].to_list()
print(mc2hpo)
hpo_dict = mc2hpo['HPO Term'].to_dict()
# Variants
snvs = pd.read_csv('../1_HPO_annotations/Annotated_files/snvs_hpo_anno.csv')
cnvs = pd.read_csv('../1_HPO_annotations/Annotated_files/cnvs_hpo_anno.csv')
strs = pd.read_csv('../1_HPO_annotations/Annotated_files/strs_hpo_anno.csv')
# Make gene symbol, gene_id map
df = cnvs[['Sample', 'Gene_Symbol', 'Gene_ID', 'GeneHPO_inclusive']]
df.columns = ['Sample', 'Gene_symbol', 'Gene_id_', 'GeneHPO_inclusive']
sym_id = pd.concat([snvs[['Sample', 'Gene_symbol', 'Gene_id_', 'GeneHPO_inclusive']], df, strs[['Sample', 'Gene_symbol', 'Gene_id_', 'GeneHPO_inclusive']]])
del df
sym_id.drop_duplicates(keep = 'first', inplace = True)

# Gene to phenotype or phenotype to genes?
def g2p(proband):
    gene = input('Enter genename:').upper()
    # Can input either an ENSEMBL ID or a gene name
    if 'ENSG' in gene:
        gene_id = gene
        gene = sym_id[sym_id.Gene_id_==gene_id]['Gene_symbol'].to_list()[0]
    else:
        gene_id = sym_id[sym_id.Gene_symbol==gene]['Gene_id_'].to_list()[0]
    print('Gene is ' + gene)
    print('Gene ID is ' + gene_id)
    # Check if proband has gene
    sub_df = sym_id[(sym_id.Sample==proband) & (sym_id.Gene_id_==gene_id)]
    if sub_df.shape[0]==0:
        print(proband + ' has no mutation in ' + gene)
        return False
    # Check for explained phenotypes
    gene_terms = sub_df.GeneHPO_inclusive.to_list()
    gene_hpo = []
    for i in gene_terms:
        gene_hpo += i.split(';')
    gene_hpo = list(set(gene_hpo))

    pro_terms = pheno[pheno.Sample==proband]['HPO_ID'].to_list()[0].split(';')

    overlap_terms = [i for i in pro_terms if i in gene_hpo]

    overlap_phenos = [hpo_dict[i] for i in overlap_terms]

    print('Explained HPO IDs are: ' + ', '.join(overlap_terms))
    print('Explained phenotypes are: ' + ', '.join(overlap_phenos))

def check_fuzzy(df):
    guess=df.iloc[0, 0]
    print("Best match is:", guess)
    right=input("Is this the correct phenotype? (Y/N)")
    if right.upper()=="Y":
        return guess
    else:
        print('Here are the next best matches:')
        print(df.iloc[1:6])
        any_right=input("Are any of these correct? (Y/N)")
        if any_right.upper()=="Y":
            index=input("What is the index of the correct phenotype? (1-5)")
            return df.iloc[str(index), 'phenotype']
        print('Please try another keyword')
        return False 

def p2g(proband):
    if proband not in pheno.Sample.to_list():
        print('Phenotype data not available for proband')
        return False
        
    phenotype=input("Enter phenotype:")
    # Can input an HPO ID or a written phenotype
    if 'HP:' in phenotype:
        hpo_term=phenotype
        phenotype=hpo_dict[hpo_term]
    else:
        # Get closest phenotype by fuzzy keyword matching
        score_info=[]
        for p in mc2hpo['MC Phenotype'].to_list():
            score=fuzz.ratio(p, phenotype)
            score_info.append([p, score])
        score_df=pd.DataFrame(score_info, columns=['phenotype', 'score'])
        score_df.sort_values('score', inplace=True, ascending=False)
        phenotype = check_fuzzy(score_df)
        if not phenotype:
            return False
        hpo_term=mc2hpo[mc2hpo['MC Phenotype']==phenotype]['HPO ID'].to_list()[0]
    print('Phenotype is:', phenotype)
    print('HPO ID is:', hpo_term)
    
    # Check if proband has phenotype
    match_lines=pheno[(pheno.Sample==proband) & (pheno.HPO_ID.str.contains(hpo_term))].shape[0]
    if match_lines==0:
        print(proband, 'does not have', phenotype)
        return False
    # Check for explained phenotypes
    genes =list(sym_id[(sym_id.Sample==proband) & (sym_id.GeneHPO_inclusive.str.contains(hpo_term))]['Gene_symbol'].unique())
    print('Second hit genes associated with phenotype:', ' '.join(genes))

while True:
    # Proband
    proband = input('Enter proband:').upper()
    print('Proband is ' + proband)
    
    selection = input('Select one of:\n(1) Enter gene and get associated phenotypes\n(2) Enter phenotype and get associated genes\n')
    if '1' in selection:
        g2p(proband)
    elif '2' in selection:
        p2g(proband)
    else:
        print('Please try again')
    
    quit=input('Do you want to check another proband? (Y/N)')
    if quit.upper()=='N':
        break