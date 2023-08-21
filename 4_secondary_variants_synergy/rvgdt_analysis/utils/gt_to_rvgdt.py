#!python3

import io
import sys
import pandas as pd

description = """

    convert genotype information to the format for RV-GDT

"""

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
        return pd.read_csv(
                io.StringIO(''.join(lines)),
                dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                       'QUAL': str, 'FILTER': str, 'INFO': str},
                sep='\t').rename(columns={'#CHROM': 'CHROM'})

def main():
    # load arguments
    vcf = sys.argv[1]
    output = sys.argv[2]

    gt = read_vcf(vcf)
    dic = { '0/0' : 0, 
            '1/0' : 1,
            '0/1' : 1,
            '1/1' : 2}
    gt = gt.replace(dic)
    cols_sample = ['SG' in i for i in list(gt.columns)]
    samples = gt.loc[:, cols_sample]
    samples_arr = samples.to_numpy().T
    with open(output, 'w') as f:
        for index, arr in enumerate(samples_arr):
            f.write(samples.columns[index] + ' ' + ' '.join([str(i) for i in arr]) + '\n')

if __name__ == '__main__':
    main()
