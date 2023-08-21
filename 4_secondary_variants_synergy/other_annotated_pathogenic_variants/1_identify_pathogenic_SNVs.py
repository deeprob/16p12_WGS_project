#!bin/python
import pandas as pd
import numpy as np

# UPDATE 7/18/2022: Definitions of pathogenicity have changed:
# 1. We will only consider LOF variants for the gene lists (missense is still ok for ClinVar)
# 2. We will filter ClinVar for the criteria provided

dropbox = 'C:/Users/corny/Dropbox/'

# This script is to identify pathogenic SNVs
snvs = pd.read_csv(dropbox+"16p12.2 project/Human patients project/WGS paper/9_Rare variant calls/Rare coding SNV calls/Rare_Deleterious_Exonic_Variants_wInheritance.txt", sep = '\t')

# Pathogenic SNVS are:
patho_snvs = snvs.copy()
merge_cols = patho_snvs.columns.to_list()
patho_snvs['CLINVAR_Pathogenic'] = '.'
# Any SNV labelled as "Pathogenic" or "likely pathogenic" in ClinVar
# AND criteria is "criteria_provided,_multiple_submitters,_no_conflicts" or "reviewed_by_expert_panel"
patho_snvs.loc[(patho_snvs.ClinVar_CLNSIG.isin(['Pathogenic', 'Likely_pathogenic', 'Pathogenic/Likely_pathogenic'])) & (patho_snvs.ClinVar_CLNREVSTAT.isin(['criteria_provided,_multiple_submitters,_no_conflicts', 'reviewed_by_expert_panel'])),
 'CLINVAR_Pathogenic'] = 'X'

# In a gene with a dominant or X-linked OMIM phenotype
omim_snv = snvs[~snvs.MIM_number.isnull()].copy()
# Get MIM translation
omim = pd.read_csv(dropbox+"16p12.2 project/Human patients project/WGS paper/5_SNV pipelines-annotations/OMIM Annotations/genemap2.txt", sep = '\t', comment = '#', header = None,
                    names = ['Chromosome', 'Genomic Position Start', 'Genomic Position End', 'Cyto Location', 'Computed Cyto Location', 'MIM Number', 'Gene Symbols', 'Gene Name',
                    'Approved Gene Symbol', 'Entrez Gene ID', 'Ensembl Gene ID', 'Comments', 'Phenotypes', 'Mouse Gene Symbol/ID'])
omim.set_index('MIM Number', inplace = True)
omim_dict = omim.Phenotypes.to_dict()
def omim_phenotype(mim):
    # If variant has only 1 mim number
    if ';' not in mim:
        return omim_dict[int(mim)]
    # If variant has multiple mim numbers
    mims = mim.split(';')
    phenos = []
    for num in mims:
        pheno = omim_dict[int(num)]
        # Skip NA
        if pheno==pheno:
            phenos.append(pheno)
    if phenos==[]:
        return np.nan
    return ';'.join(phenos)
omim_snv['OMIM_phenotype'] = omim_snv.MIM_number.apply(omim_phenotype)
print(omim_snv)
# Restrict to only autosomal dominant or X-linked genes
omim_snv = omim_snv[(omim_snv.OMIM_phenotype.str.contains('Autosomal dominant')) | (omim_snv.OMIM_phenotype.str.contains('X-linked'))]
# X-linked is only for males, so read in cohort information to get sex
cohort_info = pd.read_csv(dropbox+"16p12.2 project/Human patients project/WGS paper/3_Cohort information/16p12_All_Participants_v9.csv")
estonia = pd.read_csv(dropbox+"16p12.2 project/Human patients project/WGS paper/3_Cohort information/Estonian_Participants.csv")
# Further restrict phenotypes to remove protective and susceptability
def filter_pheno(pheno, sample):
    if sample in cohort_info.Sample.to_list():
        sex = cohort_info[cohort_info.Sample==sample]['Sex'].to_list()[0]
    else:
        sex = estonia[estonia.Sample==sample]['Sex'].to_list()[0]
    phenos = pheno.split(';')
    for dis in phenos:
        if 'protection' in dis:
            continue
        if 'susceptibility' in dis:
            continue
        if 'Autosomal dominant' in dis:
            return True
        if 'X-linked' in dis and sex=='M':
            return True
    return False
omim_snv = omim_snv[omim_snv[['OMIM_phenotype', 'Sample']].apply(lambda row: filter_pheno(row[0], row[1]), axis = 1)]
omim_snv['Dominant_X_OMIM'] = 'X'
print(omim_snv)
patho_snvs = pd.merge(patho_snvs, omim_snv, on = merge_cols, how = 'left')

# In a Tier S or 1 SFARI gene
gene_anno = pd.read_csv(dropbox+"16p12.2 project/Human patients project/WGS paper/5_SNV pipelines-annotations/Gene_Annotations/Gene_annotations.csv")
print(gene_anno.SFARI_gene_score.value_counts())
sfari_genes = gene_anno[gene_anno.SFARI_gene_score.isin(['S', '1'])][['gene_id', 'SFARI_gene_score']]
sfari_snv = pd.merge(snvs, sfari_genes, left_on = 'Gene_id_', right_on = 'gene_id', how = 'inner')
sfari_snv['SFARI_Tier_1_S'] = 'X'
print(sfari_snv)
patho_snvs = pd.merge(patho_snvs, sfari_snv, on = merge_cols, how = 'left')

# In Tier 1 or Tier 2 genes from the Geisinger DBD database
print(gene_anno.Geisinger_DBD_Tier.value_counts())
dbd_genes = gene_anno[gene_anno.Geisinger_DBD_Tier.isin(['1', '2', 1, 2])][['gene_id', 'Geisinger_DBD_Tier']]
dbd_snv = pd.merge(snvs, dbd_genes, left_on = 'Gene_id_', right_on = 'gene_id', how = 'inner')
dbd_snv['DBD_Tier_1_2'] = 'X'
print(dbd_snv)
patho_snvs = pd.merge(patho_snvs, dbd_snv, on = merge_cols, how = 'left')
print(patho_snvs[patho_snvs.Geisinger_DBD_Tier.isin(['1', '2', 1, 2])])

# Fill NA with .
patho_snvs.fillna('.', inplace = True)
# Remove extra columns
print(patho_snvs.columns)
patho_snvs = patho_snvs[['Sample', 'Chrom', 'Pos', 'Ref', 'Alt', 'Qual', 'Mut_type',
       'variant_filter', 'Func.wgEncodeGencodeBasicV19', 'Gene.wgEncodeGencodeBasicV19',
       'GeneDetail.wgEncodeGencodeBasicV19',
       'ExonicFunc.wgEncodeGencodeBasicV19',
       'AAChange.wgEncodeGencodeBasicV19', 'gnomad_exome_AF',
       'gnomad_genome_AF', 'CADD_PHRED', 'CADD_RawScore', 'ClinVar_CLNDN',
       'ClinVar_CLNDISDB', 'ClinVar_CLNREVSTAT', 'ClinVar_CLNSIG',
       'ClinVar_ALLELEID', 'GT', 'DP', 'AD', 'GQ', 'PL', 'RGQ', 'SB',
       'variant_id', 'cohort_count', 'Gene_symbol', 'Gene_id', 'Gene_biotype',
       'Gene_id_', 'LOEUF', 'MIM_number', 'inheritance', 'CLINVAR_Pathogenic',
       'OMIM_phenotype', 'Dominant_X_OMIM', 'SFARI_gene_score',
       'SFARI_Tier_1_S', 'Geisinger_DBD_Tier', 'DBD_Tier_1_2']]
print(patho_snvs)

# Add final annotation
def patho_anno(row):
    if row['CLINVAR_Pathogenic'] == 'X' or row['Dominant_X_OMIM'] == 'X' or row['SFARI_Tier_1_S'] == 'X' or row['DBD_Tier_1_2'] == 'X':
        return 'X'
    return '.'
patho_snvs['Pathogenic'] = patho_snvs.apply(lambda row: patho_anno(row), axis = 1)
print(patho_snvs.Pathogenic.value_counts())

# If variant is not LOF and not in ClinVar, remove pathogenic annotation
# Now that ClinVar are annotated, remove any non-ClinVar variant that is LOF
patho_snvs.loc[~((patho_snvs.CLINVAR_Pathogenic=='X') | (patho_snvs.Mut_type=='lof')), ['Dominant_X_OMIM', 'SFARI_Tier_1_S', 'DBD_Tier_1_2', 'Pathogenic']] = '.'

patho_snvs.to_csv('Intermediate_files/Rare_Deleterious_Exonic_Variants_Pathogenic_Anno.csv', index = False)