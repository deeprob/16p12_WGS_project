#!/bin/python3

import pandas as pd

infile = 'vcfs/protein_coding/all_chromosomes.select_annotations.tsv'
annofile = '/data5/SPARK_WES/annotations/release_june2020/gene_annotations.txt'
outfile='vcfs/protein_coding/all_chromosomes.select_annotations.leouf.tsv'

df = pd.read_csv(infile, sep='\t')

anno = pd.read_csv(annofile, sep='\t')
anno = anno.set_index('#HGNCsymbol', drop=False)
anno = anno.drop_duplicates('#HGNCsymbol')


cols = ['LOEUF', 'LOEUF_%', 'OMIM_phenotype']

genes = list(anno.index)
def get_annotation(gene, col):
        if gene in genes:
                return anno.loc[gene, col]
        return '.'

for col in cols:
        df[col] = df['Gene.refGene'].apply(lambda s: get_annotation(s, col))

df.to_csv(outfile, sep='\t', index=False)


