#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=delly
#SBATCH -o logs/3_delly_merge.log
#SBATCH -e logs/3_delly_merge.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=5G
#SBATCH --chdir /data5/16p12_WGS/structural_variants/vcf_callers_2022_02_07/delly


# bcftools merge -m id -O b -o merged.bcf s1.geno.bcf s2.geno.bcf ... sN.geno.bcf
bcftools merge -m id -O b -o /data5/16p12_WGS/structural_variants/vcf_callers_2022_02_07/delly/vcfs/3_delly_merge/merged.bcf /data5/16p12_WGS/structural_variants/vcf_callers_2022_02_07/delly/vcfs/2_regenotype/*.bcf
bcftools index /data5/16p12_WGS/structural_variants/vcf_callers_2022_02_07/delly/vcfs/3_delly_merge/merged.bcf




