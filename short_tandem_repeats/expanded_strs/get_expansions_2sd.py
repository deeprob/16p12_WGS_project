#!/bin/python3


import pandas as pd
import numpy as np
import sys

vcf = pd.read_csv('output/16p12_cohort.tsv', sep='\t', comment='#')
vcf = vcf.set_index(['chrom', 'pos'])


stats = pd.read_csv('output/16p12_cohort.statstr.tab', sep='\t')
stats = stats.set_index(['chrom', 'start'], drop=False)


bed = pd.read_csv('../vcf_callers/hg19_ver13_1.bed', sep='\t', header=None)
bed = bed.set_index([0,1], drop=False)
bed.columns = ['chrom', 'start', 'end', 'period', 'motif']
bed['ref_length'] = bed.end - bed.start + 1
def construct_ref_allele(motif, ref_length):
    return motif* int(ref_length/len(motif))

bed['ref_allele'] = [construct_ref_allele(s[0], s[1]) for s in zip(bed['motif'], bed['ref_length'])]

# remove num called == 0
# there is one entry that shouldn't be there
stats = stats[stats['numcalled'] != 0].copy()

# remove var == 0
stats = stats[stats['var'] != 0].copy()
vcf = vcf.loc[stats.index].copy()

samples = list(vcf.columns)

# write header
outline = '{}\t{}\t{}\t{}\t{}\t'.format('chrom', 'pos', 'end', 'ref_allele', 'alt_allele')
outline = outline + '{}\t{}\t{}\t{}\t{}\n'.format('sample', 'zscore', 'longest_allele', 'cohort_mode', 'motif_period')
sys.stdout.write(outline)


for i, vcf_row in vcf.iterrows():
    chrom = i[0]
    pos   = i[1]

    # get mean and std dev at locus
    loc_alleles = []
    vcf_row = vcf_row.to_dict()
    for samp in samples:
        alleles = vcf_row[samp]
        if alleles == '.':
            continue
        alleles = alleles.split(',')
        for a in alleles:
            a = int(a)
            loc_alleles.append(a)
            
    standard_dev = np.std(loc_alleles)
    mean = np.mean(loc_alleles)

    # get other stats to write out
    mode = int(stats.loc[i, 'mode'])
    end = int(stats.loc[i, 'end']) - 1
    ref_allele = bed.loc[i, 'ref_allele']
    period = bed.loc[i, 'period']
    motif = bed.loc[i, 'motif']
    
    for samp in samples:
        alleles = vcf_row[samp]
        if alleles == '.':
            continue
        alleles = alleles.split(',')
        alleles = [int(s) for s in alleles]
        longest_allele = max(alleles)
        if (longest_allele - mean) > (2 * standard_dev):
            zscore = (longest_allele - mean) / standard_dev
            alt_allele = construct_ref_allele(motif, longest_allele * period)
            outline = '{}\t{}\t{}\t{}\t{}\t'.format(chrom, pos, end, ref_allele, alt_allele)
            outline = outline + '{}\t{}\t{}\t{}\t{}\n'.format(samp, zscore, longest_allele, mode, period)
            sys.stdout.write(outline)







