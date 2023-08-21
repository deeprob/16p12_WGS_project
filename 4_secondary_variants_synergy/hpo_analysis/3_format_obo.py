#!/bin/python
import pandas as pd

# Make files with HPO relationships
obo = open('Analysis_files/hp.obo', 'r')

# Reorganize ontology
# Ontology loks like this:

# [Term]
# id: HP:0000006
# name: Autosomal dominant inheritance
# alt_id: HP:0001415
# alt_id: HP:0001447
# alt_id: HP:0001448
# alt_id: HP:0001451
# alt_id: HP:0001455
# alt_id: HP:0001456
# alt_id: HP:0001463
# def: "A mode of inheritance that is observed for traits related to a gene encoded on one of the autosomes (i.e., the human chromosomes 1-22) in which a trait manifests in heterozygotes. In the context of medical genetics, an autosomal dominant disorder is caused when a single copy of the mutant allele is present. Males and females are affected equally, and can both transmit the disorder with a risk of 50% for each child of inheriting the mutant allele." [HPO:curators]
# synonym: "Autosomal dominant" EXACT []
# synonym: "Autosomal dominant form" RELATED [HPO:skoehler]
# synonym: "Autosomal dominant type" RELATED [HPO:skoehler]
# xref: SNOMEDCT_US:263681008
# xref: UMLS:C0443147
# is_a: HP:0000005 ! Mode of inheritance

# For each entry in the file, save all information to a dataframe
newrows = []
rel_annos = ['id', 'name', 'alt_id', 'def', 'synonym', 'comment', 'consider', 'created_by', 'creation_date', 'is_a', 'is_obsolete', 'property_value', 'replaced_by', 'subset', 'xref']
anno_dict = {}
for i in rel_annos:
    anno_dict[i] = []
    
for line in obo:
    line = line.rstrip()
    if line=='[Term]':
        if anno_dict['id']!=[]:
            newrows.append([';'.join(anno_dict[i]) for i in rel_annos])
        for i in rel_annos:
            anno_dict[i] = []
        
        continue
    
    line = line.split(': ')
    for anno in rel_annos:
        if line[0]==anno:
            anno_dict[anno].append(line[1])
            continue
newrows.append([';'.join(anno_dict[i]) for i in rel_annos])

# Convert to DataFrame
df = pd.DataFrame(newrows, columns = rel_annos)
df.replace('', '.', inplace = True)
print(df)

# Save to file
df.to_csv('Analysis_files/3_hpo_obo_reformat.csv', index = False)

