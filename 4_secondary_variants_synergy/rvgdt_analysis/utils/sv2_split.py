#!/usr/bin/env python

description="Split the structural variant calls to individual files"

import sys
import pandas as pd
from multiprocessing import Pool

intermediate_files_dir="/data5/deepro/wgs_16p/rvgdt/data/intermediate_files"

# this function was added by Deepro since the format of the new sv calls
# file has changed
def convert_to_bed_like_format(row):
    chrm, start, end = row.call_ID.split("_")[:3]
    return pd.Series([chrm, int(start), int(end)], index=["chrom", "start", "end"])

def save(df):
    """
    for each piece of dataframe, save it by sample ID
    """
    bed = df.apply(convert_to_bed_like_format, axis=1) # edited by Deepro; new sv calls file format
    bed.start -= 1
    bed.to_csv(f"{intermediate_files_dir}/svs/{df.iloc[0,0]}.sv2.bed", header=False, index=False, sep="\t")
 
def main():
    sv = pd.read_csv(sys.argv[1], sep="\t", )
    dfs = sv.groupby('Sample')  # group dataframe by sample IDs
    dfs = [pd.DataFrame(y) for x, y in dfs] # list comprehension of dataframes
   
    with Pool(10) as p:
        p.map(save, dfs)
   
if __name__ == '__main__':
    main()
