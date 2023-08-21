#!/bin/python
import pandas as pd

dropbox = 'D:/Dropbox/'

# Annotate variant files with HPO terms for easy access

# SNVs
snvs = pd.read_csv(dropbox+"16p12.2 project/Human patients project/WGS paper/12_Pathogenic Variant Analysis/1_Pathogenic_SNV/Rare_Deleterious_Exonic_Variants_Pathogenic_Anno.csv")
print(snvs)

# Annotate with HPO terms
gene2hpo = pd.read_csv('Analysis_files/2_ensembl2hpo.csv')
gene2hpo.set_index(gene2hpo.ENSEMBL_ID, inplace = True)
gene2hpo_dict = gene2hpo.HPO_IDs.to_dict()

snvs['Gene2HPO'] = snvs.Gene_id_.map(gene2hpo_dict)
snvs.fillna('.', inplace = True)

# Annotate with alternate terms of gene HPO terms
hpo_onto = pd.read_csv('Analysis_files/5_term_relationships.csv')

def anno_alternate(hpo):
    if hpo=='.':
        return '.'
    ids = hpo.split(';')
    newterms = hpo_onto[hpo_onto.id.isin(ids)]['alt_id'].to_list()
    outterms = []
    for i in newterms:
        if i!='.':
            outterms += i.split(';')
    output = list(set(ids + outterms))
    output.sort()
    return ';'.join(output)

snvs['GeneHPO_inclusive'] = snvs.Gene2HPO.apply(anno_alternate)

def anno_parent(hpo, parent = 'direct'):
    if hpo=='.':
        return '.'
    ids = hpo.split(';')
    if parent == 'direct':
        newterms = hpo_onto[hpo_onto.id.isin(ids)]['direct_parent_ids'].to_list()
    elif parent == 'indirect':
        newterms = hpo_onto[hpo_onto.id.isin(ids)]['indirect_parent_ids'].to_list()
    outterms = []
    for i in newterms:
        if i!='.':
            outterms += i.split(';')
    output = list(set(outterms))
    output.sort()
    return ';'.join(output)

snvs['HPO_parent_direct'] = snvs.GeneHPO_inclusive.apply(anno_parent)
snvs['HPO_parent_indirect'] = snvs.GeneHPO_inclusive.apply(anno_parent, parent = 'indirect')
print(snvs)

# Save to file
snvs.to_csv('Annotated_files/snvs_hpo_anno.csv', index = False)

# CNVs
cnvs = pd.read_csv(dropbox+"16p12.2 project/Human patients project/WGS paper/9_Rare variant calls/Rare coding CNV calls/sv_calls_combined.txt", sep = '\t')
cnvs['Gene2HPO'] = cnvs.Gene_ID.map(gene2hpo_dict)
cnvs.fillna('.', inplace = True)
cnvs['GeneHPO_inclusive'] = cnvs.Gene2HPO.apply(anno_alternate)
cnvs['HPO_parent_direct'] = cnvs.GeneHPO_inclusive.apply(anno_parent)
cnvs['HPO_parent_indirect'] = cnvs.GeneHPO_inclusive.apply(anno_parent, parent = 'indirect')
print(cnvs)
cnvs.to_csv('Annotated_files/cnvs_hpo_anno.csv', index = False)

strs = pd.read_csv(dropbox+"16p12.2 project/Human patients project/WGS paper/8_STR variants/Exonic_STRs_wInheritance.csv")
strs['Gene2HPO'] = strs.Gene_id_.map(gene2hpo_dict)
strs.fillna('.', inplace = True)
strs['GeneHPO_inclusive'] = strs.Gene2HPO.apply(anno_alternate)
strs['HPO_parent_direct'] = strs.GeneHPO_inclusive.apply(anno_parent)
strs['HPO_parent_indirect'] = strs.GeneHPO_inclusive.apply(anno_parent, parent = 'indirect')
print(strs)
strs.to_csv('Annotated_files/strs_hpo_anno.csv', index = False)





