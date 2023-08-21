#!/bin/python
import pandas as pd

# Assess enrichment of genes in different BrainSpan timelines
df = pd.read_csv('Analysis_files/brainspan_preferential_tissue_expression_minor_epoch.tsv', sep = '\t')

# Get list of all genes that both could be expressed and could be variant in our cohort
exp_genes = list(df.ensembl_gene_id.unique())
var_genes = list(pd.read_csv('../../5_SNV pipelines-annotations/Gene_Annotations/Gene_annotations.csv')['gene_id'].unique())
all_genes = list(set(exp_genes) & set(var_genes))


# Get gene lists for our cohort
snvs = pd.read_csv('../../9_Rare variant calls/Rare coding SNV calls/Rare_Deleterious_Exonic_Variants.csv')
cnvs = pd.read_csv('../../9_Rare variant calls/Rare coding CNV calls/sv_calls_combined.txt', sep = '\t')
cnvs['Mut_type'] = cnvs.Type
strs = pd.read_csv('../../8_STR variants/Exonic_STR_expansions.csv')
strs.Mut_type = 'STR'

# Combine all variants into one file
sub_snv = snvs[['Sample', 'Mut_type', 'Gene_id_', 'Gene_symbol', 'variant_id']]
sub_cnv = cnvs[['Sample', 'Type', 'Gene_ID', 'Gene_Symbol', 'call_ID']]
sub_cnv.columns = ['Sample', 'Mut_type', 'Gene_id_', 'Gene_symbol', 'variant_id']
sub_str = strs[['Sample', 'Gene_id_', 'Gene_symbol', 'variant_id']].copy()
sub_str['Mut_type'] = 'STR'
var_df = pd.concat([sub_snv, sub_cnv, sub_str], ignore_index = True)

# Save memory
del snvs
del cnvs
del strs
del sub_snv
del sub_cnv
del sub_str

# Restrict to only genes in probands or their parents
cohort_info = pd.read_csv('../../3_Cohort information/16p12_All_Participants_v9.csv')
cohort_info = cohort_info[(cohort_info.WGS=='X') & (cohort_info.No_consent_forms!='X')]

relation_dict = dict(zip(cohort_info.Sample.to_list(), cohort_info.Relationship.to_list()))
var_df['Relationship'] = var_df.Sample.map(relation_dict)

var_df = var_df[var_df.Relationship.isin(['P', 'MC', 'MNC', 'FC', 'FNC'])]

print(var_df)

# Get numbers for enrichment calculations
# a) Expressed and has variant in our cohort
# b) Expressed and does not have variant in our cohort
# c) Not expressed and has variant in our cohort
# d) Not expressed and does not have variant in our cohort

# Need to calculate these values per region per time per variant type per relationship status
regions = ['cerebellum', 'hippocampus', 'amygdala', 'striatum', 'thalamus', 'medial frontal cortex', 'orbitofrontal cortex',
           'dorsolateral frontal cortex', 'ventrolateral frontal cortex', 'primary motor cortex',
           'primary somatosensory cortex', 'inferior parietal cortex', 'primary auditory cortex',
           'superior temporal cortex', 'inferior temporal cortex', 'primary visual cortex']
times = ['early fetal', 'early mid-fetal', 'late mid-fetal', 'late fetal', 'early infancy', 'late infancy',
         'early childhood', 'middle and late childhood', 'adolescence', 'young adulthood', 'middle adulthood']
count_calcs = []
rel_code = {'proband':['P'], 'carrier_parent':['MC', 'FC'], 'noncarrier_parent':['MNC', 'FNC']}
for rel in ['proband', 'carrier_parent', 'noncarrier_parent']:
    for region in regions:
        for time in times:
            # Make sure region-time combination exists
            rt_cols = df.columns[(df.columns.str.contains(region)) & (df.columns.str.contains(time))].to_list()
            if len(rt_cols)==0:
                print(region, time)
                continue
            # Get list of genes expressed at region-time
            rt_genes = list(set(all_genes) & set(df[df[rt_cols[0]]==1]['ensembl_gene_id'].to_list()))
            for type in ['missense', 'lof', 'splice', 'DEL', 'DUP', 'STR']:
                # Get genes of type in relationship
                type_genes = list(set(all_genes) & set(var_df[(var_df.Relationship.isin(rel_code[rel])) & (var_df.Mut_type==type)]['Gene_id_'].to_list()))
                # Get counts of abcd (above) for type, region, time, relationship
                expr_var = len(list(set(rt_genes) & set(type_genes)))
                expr_novar = len(list(set(rt_genes) - set(type_genes)))
                noexpr_var = len(list(set(type_genes) - set(rt_genes)))
                noexpr_novar = len(list((set(all_genes) - set(type_genes)) - set(rt_genes)))
                
                output = [rel, region, time, type, expr_var, expr_novar, noexpr_var, noexpr_novar]
                count_calcs.append(output)
            
            # Also look at all genes, regardless of variant type
            type_genes = var_df[var_df.Relationship.isin(rel_code[rel])]['Gene_id_'].to_list()
            
            expr_var = len(list(set(rt_genes) & set(type_genes)))
            expr_novar = len(list(set(rt_genes) - set(type_genes)))
            noexpr_var = len(list(set(type_genes) - set(rt_genes)))
            noexpr_novar = len(list((set(all_genes) - set(type_genes)) - set(rt_genes)))
                
            output = [rel, region, time, 'All', expr_var, expr_novar, noexpr_var, noexpr_novar]
            count_calcs.append(output)
            print(rel, region, time)

count_df = pd.DataFrame(count_calcs, columns = ['relationship', 'region', 'time', 'type', 'expressed_variant', 'expressed_novariant', 'notexpressed_variant', 'notexpressed_novariant'])
print(count_df)
count_df.to_csv('Analysis_files/1_expression_variant_counts.csv', index = False)