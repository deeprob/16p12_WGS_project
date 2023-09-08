#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=delly
#SBATCH -o logs/5_filter.log
#SBATCH -e logs/5_filter.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=5G
#SBATCH --chdir /data5/16p12_WGS/structural_variants/vcf_callers_2022_02_07/delly



echo `date` job started on $HOSTNAME


export PATH=/data5/anastasia/sw/bcftools-1.12:$PATH

infile=/data5/16p12_WGS/structural_variants/vcf_callers_2022_02_07/delly/vcfs/4_delly_filter/merged.bcf
outfile=/data5/16p12_WGS/structural_variants/vcf_callers_2022_02_07/delly/vcfs/5_filter/merged.vcf.gz


# filter with locus-specific (FILTER==PASS) and sample-specific (FT==PASS) default filters
# when filtering for locus-specific, use 
# 	bcftools view -i 'FILTER="PASS" & FORMAT/FT="PASS"' $invcf
# This will remove records where the variant itself is low quality and where no samples have a variant that passes.
# But genotype information for Samples with FT=lowQuality will still be present
# to remove that information, use
# 	bcftools filter --set-GTs . -i 'FORMAT/FT="PASS"'
# to set genotypes to missing (GT="./.") for low quality variants


bcftools view -i 'FILTER="PASS" & FORMAT/FT="PASS" & FORMAT/GT="alt"' $infile | bcftools filter --set-GTs . -i 'FORMAT/FT="PASS" & FORMAT/GT="alt"' |  bgzip > $outfile


# inspect output with
# bcftools query -f '%FILTER\t[%GT:%FT\t]\n' $infile | less -S
# bcftools query -f '%FILTER\t[%GT:%FT\t]\n' $outfile | less -S


echo `date` job finished



