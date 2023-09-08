#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=delly_cnvs
#SBATCH -o logs/8_str_filter.log
#SBATCH -e logs/8_str_filter.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --chdir /data5/16p12_WGS/structural_variants/sv_caller_postprocessing/2022_02_14/delly/


# use bedtools to find either
# (a) CNVs that do not overlap the STR regions
# <-------->                <---------------> CNV
#             <---->                          STR
bedtools intersect -wa -a bed_files/7_brandler_filter.bed -b ../hg19_ver13_1.bed -v -header > bed_files/8.1_str_no_overlap.bed
# (b) CNVS that overlap the STR regions 100%
# <-----------------------> CNV
#    <--->                  STR
bedtools intersect -wa -a bed_files/7_brandler_filter.bed -b ../hg19_ver13_1.bed -F 1 -u -header > bed_files/8.2_str_all_overlap.bed

# Anastasia noticed a problem where this second line would also get:
# <---------------------->    CNV
#    <---->            <----> STR
# So run a python script to filter the 100% overlap calls for those with a breakpoint in an STR
python 8.1_str_filter.py bed_files/8.2_str_all_overlap.bed

# Merge two files
cat bed_files/8.1_str_no_overlap.bed > bed_files/8.4_str_brandler_filter.bed
cat bed_files/8.3_all_str_filter.bed >> bed_files/8.4_str_brandler_filter.bed

# Add a header
head -1 bed_files/8.1_str_no_overlap.bed > bed_files/8_str_brandler_filter.bed
tail -n+2 bed_files/8.4_str_brandler_filter.bed | sort | uniq >> bed_files/8_str_brandler_filter.bed

