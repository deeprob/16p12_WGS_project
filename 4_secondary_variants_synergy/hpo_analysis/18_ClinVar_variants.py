import pandas as pd

# Get ClinVar variants in probands for manual HPO annotation

# Load in data
# Get probands
cohort_info = pd.read_csv('../../3_cohort information/16p12_All_Participants_v9.csv')
cohort_info = cohort_info[(cohort_info.WGS=='X') & (cohort_info.No_consent_forms!='X')]
pros = cohort_info[cohort_info.Relationship=='P']['Sample'].to_list()

# Variants
snvs = pd.read_csv('../1_HPO_annotations/Annotated_files/snvs_hpo_anno.csv')
snvs = snvs[(snvs.Sample.isin(pros)) & (snvs.CLINVAR_Pathogenic=='X')]

# Save to file
snvs.to_csv('Analysis_files/4_ClinVar_variants.csv', index=False)