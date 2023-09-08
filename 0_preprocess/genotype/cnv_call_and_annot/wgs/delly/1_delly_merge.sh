#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=delly
#SBATCH -o logs/1_merge.log
#SBATCH -e logs/1_merge.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=5G
#SBATCH --chdir /data5/16p12_WGS/structural_variants/vcf_callers_2022_02_07/delly



echo `date` starting job on $HOSTNAME as $USER




delly=/data3/software/delly_v0.8.2_linux_x86_64bit


# delly merge -o sites.bcf s1.bcf s2.bcf ... sN.bcf
$delly merge -o /data5/16p12_WGS/structural_variants/vcf_callers_2022_02_07/delly/vcfs/1_merge/merged.bcf /data5/16p12_WGS/structural_variants/vcf_callers/tmp/delly/*.bcf




echo `date` finished




