#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=delly
#SBATCH -o logs/4_delly_filter.log
#SBATCH -e logs/4_delly_filter.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=5G
#SBATCH --chdir /data5/16p12_WGS/structural_variants/vcf_callers_2022_02_07/delly



delly=/data3/software/delly_v0.8.2_linux_x86_64bit


# delly filter -f germline -o germline.bcf merged.bcf
$delly filter -f germline -o /data5/16p12_WGS/structural_variants/vcf_callers_2022_02_07/delly/vcfs/4_delly_filter/merged.bcf /data5/16p12_WGS/structural_variants/vcf_callers_2022_02_07/delly/vcfs/3_delly_merge/merged.bcf





