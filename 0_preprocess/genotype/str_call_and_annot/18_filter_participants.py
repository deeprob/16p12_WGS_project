#!/bin/python3


# combine tables


import pandas as pd


# 16p12_All_Participants_v5.csv is from /Dropbox/16p12.2 project/Human patients project/WGS paper/3_Cohort information/16p12_All_Participants_v5.csv

samples_df = pd.read_csv('16p12_All_Participants_v5.csv')
samples_df = samples_df[samples_df['No_consent_forms'] != 'X']
samples_df = samples_df[samples_df['WGS'] == 'X']




df = pd.read_csv('16p12_cohort.expansions_2SD.annotated.exonic.gene_ids.column_names.loeuf.omim.csv')

df = df[df['Sample'].isin(samples_df['Sample'])]



df.to_csv('16p12_cohort.expansions_2SD.annotated.exonic.gene_ids.column_names.loeuf.omim.participants.csv')
















