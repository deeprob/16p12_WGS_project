#!/bin/python3


import sys
import pandas as pd

filename = sys.argv[1]
outfilename = sys.argv[2]

locs = pd.read_csv('../vcf_callers/hg19_ver13_1.bed', sep='\t', header=None, names=['chrom', 'start', 'end', 'motif_length', 'motif'])
locs['var_id'] = locs['chrom'].apply(lambda s: s[3:]) + '_' + locs['start'].astype(str)
locs = locs.set_index('var_id', drop=False)


de_novo = pd.read_csv(filename, sep='\t')

str_ids = de_novo.str_id.to_list()

# remove ids that aren't de novo
locs = locs[locs.var_id.isin(str_ids)]


locs['ref_length'] = locs.end - locs.start + 1

def construct_ref_allele(motif, ref_length):
    return motif* int(ref_length/len(motif))
locs['ref_allele'] = [construct_ref_allele(s[0], s[1]) for s in zip(locs['motif'], locs['ref_length'])]

# dummy alt allele
locs['alt_allele'] = locs['ref_allele'].apply(lambda s: s[:-1])

cols = ['chrom', 'start', 'end', 'ref_allele', 'alt_allele']

locs = locs[cols]


locs.to_csv(outfilename, sep='\t', index=False, header=False)

