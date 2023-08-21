#!python3

description = """
    Map phenotype split to the family pedigree file to do RV-GDT analysis (prelim_child_domains_Jul16.xlsx)

    $ python map_phenotype_split.py
"""
import os
import numpy as np
import pandas as pd

### Added dir and filenames ### 
INPUT_DIR="/data5/deepro/wgs_16p/rvgdt/data/input_files"
INTERMEDIATE_DIR="/data5/deepro/wgs_16p/rvgdt/data/intermediate_files"

DOMAIN_INFO_FILE=f"{INPUT_DIR}/prelim_child_domains_Jul16.xlsx"
PHENO_SPLITS_FILE=f"{INPUT_DIR}/Phenotype_splits_new.xlsx"
PED_FILE=f"{INPUT_DIR}/Mar_2022.fam"
### ------------------------ ###

if __name__ == '__main__':
    # load the phenotype score and the corresponding split values
    phenotypes = pd.read_excel(DOMAIN_INFO_FILE, header=None, engine='openpyxl') # changed to variable
    splits = pd.read_excel(PHENO_SPLITS_FILE, header=0, engine='openpyxl') # changed to variable
    
    # transpose the dataframe and set sample ID as index
    phenotypes.set_index(0, inplace=True)
    phenotypes = phenotypes.transpose()
    phenotypes.set_index("Child Code", inplace=True)

    # keep the domain columns
    phenotypes = phenotypes.iloc[:, [3, 8, 14, 20, 24, 31]]
    phenotypes.columns = ["Child_ID_DD", "Child_behav", "Child_psych", "Child_nervous_system", "Child_congenital", "Child_craniofacial"]

    # replace NA with 0, all values > 0 are replaced with 0
    splits.set_index("Phenotype", inplace=True, drop=True)
    phenotypes.fillna(0, inplace=True)
    for p in phenotypes.columns:
        phenotypes.loc[:, [p]] = np.select([phenotypes.loc[:, [p]] > splits.loc[p, "Split_value"]], [1], default=0) 

    # For each binary phenotype, map them back to the pedigree file, samples without phenotype info will be marked as 0
    family = pd.read_table(PED_FILE, header=0) # edited header to 0 also changed ped file to newer fam file
    family.columns = ["Family", "Sample", "Father", "Mother", "Sex", "Carrier"]
    family = family.loc[:, ["Family", "Sample", "Father", "Mother", "Sex"]]

    ### Adding dir creation command ###
    os.makedirs(os.path.join(INTERMEDIATE_DIR, "pedigrees"))
    ### --------------------------- ###

    for p in phenotypes.columns:
        df = family.copy()
        df = pd.merge(df, phenotypes.loc[:, [p]], left_on = "Sample", right_index = True, how = "left")
        df.fillna(0, inplace=True)
        df = df.astype({p: int})
        # unaffected 1, affected 2
        df.loc[:, [p]] = df.loc[:, [p]] + 1
        df.to_csv(f"{INTERMEDIATE_DIR}/pedigrees/family_pedigree.{p}.txt", header=None, sep="\t", index=False) # added dir prefix
