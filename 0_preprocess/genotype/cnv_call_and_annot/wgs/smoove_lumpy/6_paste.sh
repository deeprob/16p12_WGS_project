#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=smoove
#SBATCH -o logs/4_paste.log
#SBATCH -e logs/4_paste.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=5G
#SBATCH --chdir /data5/16p12_WGS/structural_variants/vcf_callers_2022_02_07/smoove_lumpy
#SBATCH --nodelist qingyu


echo `date` starting job on $HOSTNAME as $USER


fasta=/data5/bx_reference/hg19/ucsc.hg19.fasta
in_dir=/data5/16p12_WGS/structural_variants/vcf_callers_2022_02_07/smoove_lumpy/vcfs/3_regenotype/
out_dir=/data5/16p12_WGS/structural_variants/vcf_callers_2022_02_07/smoove_lumpy/vcfs/4_paste/



export PATH=/data5/anastasia/sw/bcftools-1.12:$PATH
export PATH=/data3/software/bin:$PATH
export PATH=/data5/software/lumpy-sv/bin/:$PATH

# Instruction at https://github.com/brentp/smoove

smoove paste --name merged --outdir $out_dir $in_dir




echo `date` finished




