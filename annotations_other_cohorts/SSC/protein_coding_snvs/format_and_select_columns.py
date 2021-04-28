#!/bin/python3

# this script does the final formatting and selects annotations for final table

import sys
import pandas as pd


infilename = 'vcfs/protein_coding/all_chromosomes.txt'
sample_map_filename = '../documentation/nygc_sfari_id_map.csv'
outfilename = 'vcfs/protein_coding/all_chromosomes.select_annotations.tsv'

# load in samp mapp
sampmapp = pd.read_csv(sample_map_filename)
sampmapp = sampmapp.set_index('Repository Id')['SFARI ID'].to_dict()


df = pd.read_csv(infilename, sep='\t', header=None)

# name original columns
df.columns = ['Chr', 'Pos', 'ID', 'Ref', 'Alt', 'Qual', 'Filter', 'Info', 'SampleID']

# annotations we want
cols = ['Family','Sample','Chr','Pos','Ref','Alt','Qual', 'Mut_type',
        'Gene.refGene', 'Num_genes', 'RVIS_%.refGene', 'pLI.refGene',
        'CADD', 'gnomAD', 'MPC', 'GERP++_RS', 'GERP++_NR', 'ASD_risk_genes_TADA_FDR0.3',
        'ASD_coexpression_networks_Willsey2013','BrainExpressed_Kang2011',
        'PSD_Genes2Cognition','Developmental_delay_DDD',
        'CHD8_targets_Cotney2015_Sugathan2014','FMRP_targets_Darnell2011',
        'DBD','DDDG2P','SFARI_Gene','Purcell_Schiz']

def get_info_field(info_line, info_field):
        # returns the info field
        sinfo = info_line.split(';')
        for item in sinfo:
                if item.startswith(info_field):
                        length_info_field = len(info_field)
                        return item[length_info_field + 1:]
        return '.'

# make all the columns
df['Sample'] = df.SampleID.apply(lambda s: sampmapp[s])
df['Family'] = df.Sample.apply(lambda s: s.split('.')[0])
df['Mut_type'] = df.Info.apply(lambda s: get_info_field(s, 'variant_function'))
df['Gene.refGene'] = df.Info.apply(lambda s: get_info_field(s, 'Gene.refGene'))
df['Gene.refGene'] = df['Gene.refGene'].apply(lambda s: s.split('\\x3b'))
df['Gene.refGene'] = df['Gene.refGene'].apply(lambda s: ';'.join(s))
df['Num_genes'] = df['Gene.refGene'].apply(lambda s: len(s.split(';')))
df['RVIS_%.refGene'] = df.Info.apply(lambda s: get_info_field(s, 'RVIS_%.refGene'))
df['pLI.refGene'] = df.Info.apply(lambda s: get_info_field(s, 'pLI.refGene'))
df['CADD'] = df.Info.apply(lambda s: get_info_field(s, 'CADD_phred'))
df['gnomAD'] = df.Info.apply(lambda s: get_info_field(s, 'gnomad_AF'))
df['MPC'] = df.Info.apply(lambda s: get_info_field(s, 'MPC'))
df['GERP++_RS'] = df.Info.apply(lambda s: get_info_field(s, 'GERP++_RS'))
df['GERP++_NR'] = df.Info.apply(lambda s: get_info_field(s, 'GERP++_NR'))
df['ASD_risk_genes_TADA_FDR0.3'] = df.Info.apply(lambda s: get_info_field(s, 'ASD_risk_genes_TADA_FDR0.3.refGene'))
df['ASD_coexpression_networks_Willsey2013'] = df.Info.apply(lambda s: get_info_field(s, 'ASD_coexpression_networks_Willsey2013.refGene'))
df['BrainExpressed_Kang2011'] = df.Info.apply(lambda s: get_info_field(s, 'BrainExpressed_Kang2011.refGene'))
df['PSD_Genes2Cognition'] = df.Info.apply(lambda s: get_info_field(s, 'PSD_Genes2Cognition.refGene'))
df['Developmental_delay_DDD'] = df.Info.apply(lambda s: get_info_field(s, 'Developmental_delay_DDD.refGene'))
df['CHD8_targets_Cotney2015_Sugathan2014'] = df.Info.apply(lambda s: get_info_field(s, 'CHD8_targets_Cotney2015_Sugathan2014.refGene'))
df['FMRP_targets_Darnell2011'] = df.Info.apply(lambda s: get_info_field(s, 'FMRP_targets_Darnell2011.refGene'))
df['DBD'] = df.Info.apply(lambda s: get_info_field(s, 'DBD.refGene'))
df['DDDG2P'] = df.Info.apply(lambda s: get_info_field(s, 'DDDG2P.refGene'))
df['SFARI_Gene'] = df.Info.apply(lambda s: get_info_field(s, 'SFARI_Gene.refGene'))
df['Purcell_Schiz'] = df.Info.apply(lambda s: get_info_field(s, 'Purcell_Schiz.refGene'))

# keep the columns we want
df = df[cols]

# add inheritence
df['variant_id'] = df.Chr + '_' + df.Pos.astype(str) + '_' + df.Alt
df['inheritence'] = '.'

total_num_families = len(df.Family.unique())
num_family = 0

for fam, fam_df in df.groupby('Family'):
        print('{} {}/{}'.format(fam, num_family, total_num_families))
        num_family = num_family + 1
        # if either parent is not in file then skip
        samples_in_family = list(fam_df.Sample.unique())
        if str(fam) + '.mo' not in samples_in_family:
                continue
        if str(fam) + '.fa' not in samples_in_family:
                continue
        father_variants = fam_df[fam_df.Sample == str(fam) + '.fa'].variant_id.to_list()
        mother_variants = fam_df[fam_df.Sample == str(fam) + '.mo'].variant_id.to_list()
        for i, row in fam_df.iterrows():
                samp = row['Sample']
                # if it's one of the parents then skip
                if samp == (str(fam) + '.mo') or samp == (str(fam) + '.fa'):
                        continue
                var_id = row['variant_id']
                if var_id in father_variants and var_id in mother_variants:
                        inheritence = 'both'
                elif var_id in father_variants:
                        inheritence = 'father'
                elif var_id in mother_variants:
                        inheritence = 'mother'
                else:
                        inheritence = 'de novo'
                df.at[i, 'inheritence'] = inheritence



# save to file
df.to_csv(outfilename, sep='\t', index=False)




