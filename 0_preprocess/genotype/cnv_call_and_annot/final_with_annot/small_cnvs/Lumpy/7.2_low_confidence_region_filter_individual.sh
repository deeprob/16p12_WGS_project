#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=low_confidence_filter
#SBATCH -o logs/7_brandler_filter.log
#SBATCH -e logs/7_brandler_filter.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=2G
#SBATCH --chdir /data5/16p12_WGS/structural_variants/sv_caller_postprocessing/2022_03_02/lumpy

# use bedtools to find regions that have less than a 50% reciprocal overlap with centromeres, SegDups, regions of low mappability, and V(D)J recombination regions from Brandler et al. Science 2016
for i in {2..346};
do
	SAMPLE=`head -n $i lumpy_files.csv | tail -n 1 | cut -f 2 -d ,`
	echo $SAMPLE `date`
	bedtools intersect -wa -a bed_files/4_size_filter/$SAMPLE'_size_filter.bed' -b ../../2022_02_14/filtered_regions_hg19.bed -f 0.5 -r -v -header > bed_files/7_brandler_filter/$SAMPLE'_brandler_filter.bed'
done
