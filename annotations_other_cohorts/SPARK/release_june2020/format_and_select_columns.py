#!/bin/python3

# this script does the final formatting and selects annotations for final table

import sys
import pandas as pd


infilename = 'vcfs/formatted/ssc_release_june_2020.txt'
sample_map_filename = '/data5/SPARK_WES/Mastertable/SPARK.WES2.mastertable.2021_03.tsv'
outfilename = 'vcfs/formatted/ssc_release_june_2020.select_annotations.txt'

# load in samp mapp
sampmapp = pd.read_csv(sample_map_filename, sep='\t')
sampmapp = sampmapp.set_index('spid')


df = pd.read_csv(infilename, sep='\t', header=None)

# name original columns
df.columns = ['Chr', 'Pos', 'ID', 'Ref', 'Alt', 'Qual', 'Filter', 'Info', 'Format', 'gt_info', 'Sample']

# annotations we want
cols = ['Family','Sample','Chr','Pos','Ref','Alt','Qual', 'Func.refGene','ExonicFunc.refGene',
        'Gene.refGene', 'Num_genes', 'RVIS_%.refGene', 'pLI.refGene',
        'CADD', 'gnomAD', 'MPC', 'ASD_risk_genes_TADA_FDR0.3',
        'ASD_coexpression_networks_Willsey2013','BrainExpressed_Kang2011',
        'PSD_Genes2Cognition','Developmental_delay_DDD',
        'CHD8_targets_Cotney2015_Sugathan2014','FMRP_targets_Darnell2011',
        'DBD','DDDG2P','SFARI_Gene','Purcell_Schiz', 'LOEUF', 'OMIM_phenotype']

def get_info_field(info_line, info_field):
        # returns the info field
        sinfo = info_line.split(';')
        for item in sinfo:
                if item.startswith(info_field + '='):
                        length_info_field = len(info_field)
                        return item[length_info_field + 1:]
        return '.'


# make all the columns
df['Family'] = df.Sample.apply(lambda s: sampmapp.loc[s, 'sfid'])
df['Func.refGene'] = df.Info.apply(lambda s: get_info_field(s, 'Func.refGene'))
df['ExonicFunc.refGene'] = df.Info.apply(lambda s: get_info_field(s, 'ExonicFunc.refGene'))
df['Gene.refGene'] = df.Info.apply(lambda s: get_info_field(s, 'Gene.refGene'))
df['Gene.refGene'] = df['Gene.refGene'].apply(lambda s: s.split('\\x3b'))
df['Gene.refGene'] = df['Gene.refGene'].apply(lambda s: ';'.join(s))
df['Num_genes'] = df['Gene.refGene'].apply(lambda s: len(s.split(';')))
df['RVIS_%.refGene'] = df.Info.apply(lambda s: get_info_field(s, 'RVIS_%.refGene'))
df['pLI.refGene'] = df.Info.apply(lambda s: get_info_field(s, 'pLI.refGene'))
df['CADD'] = df.Info.apply(lambda s: get_info_field(s, 'PHRED'))
df['gnomAD'] = df.Info.apply(lambda s: get_info_field(s, 'gnomad_AF'))
df['MPC'] = df.Info.apply(lambda s: get_info_field(s, 'MPC'))
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
df['LOEUF'] = df.Info.apply(lambda s: get_info_field(s, 'LOEUF.refGene'))
df['OMIM_phenotype'] = df.Info.apply(lambda s: get_info_field(s, 'OMIM_phenotype.refGene'))



# keep the columns we want
df = df[cols]

# add inheritence
df['variant_id'] = df.Chr.astype(str) + '_' + df.Pos.astype(str) + '_' + df.Alt
df['inheritence'] = '.'

fathers = list(sampmapp['father'].unique())
mothers = list(sampmapp['mother'].unique())


total_num_families = len(df.Family.unique())
num_family = 0


for fam, fam_df in df.groupby('Family'):
    print('{} {}/{}'.format(fam, num_family, total_num_families))
    num_family = num_family + 1
    #
    samples_in_family = list(fam_df.Sample.unique())
    #
    # find out who is father and who is mother
    mother = '.'
    father = '.'
    for samp in samples_in_family:
        if samp in mothers:
            mother = samp
        if samp in fathers:
            father = samp
    #
    # if either mother or father not in file then skip
    if mother == '.' or father == '.':
        continue
    # 
    # get variants of father and mother
    father_variants = fam_df[fam_df.Sample == father].variant_id.to_list()
    mother_variants = fam_df[fam_df.Sample == mother].variant_id.to_list()
    #
    # loop through each sample
    for i, row in fam_df.iterrows():
        samp = row['Sample']
        #
        # if it's one of the parents then skip
        if samp == father or samp == mother:
                continue
        #
        var_id = row['variant_id']
        #
        if var_id in father_variants and var_id in mother_variants:
                inheritence = 'both'
        elif var_id in father_variants:
                inheritence = 'father'
        elif var_id in mother_variants:
                inheritence = 'mother'
        else:
                inheritence = 'de novo'
        #
        df.at[i, 'inheritence'] = inheritence






# save to file
df.to_csv(outfilename, sep='\t', index=False)




