#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=smoove
#SBATCH -o logs/2_union.log
#SBATCH -e logs/2_union.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=5G
#SBATCH --chdir /data5/16p12_WGS/structural_variants/vcf_callers_2022_02_07/smoove_lumpy
#SBATCH --nodelist qingyu



echo `date` starting job on $HOSTNAME as $USER


fasta=/data5/bx_reference/hg19/ucsc.hg19.fasta
in_dir=/data5/16p12_WGS/structural_variants/vcf_callers_2022_02_07/smoove_lumpy/vcfs/1_smoove/
out_dir=/data5/16p12_WGS/structural_variants/vcf_callers_2022_02_07/smoove_lumpy/vcfs/2_union/



export PATH=/data5/anastasia/sw/bcftools-1.12:$PATH
export PATH=/data3/software/bin:$PATH
export PATH=/data5/software/lumpy-sv/bin/:$PATH


# Instruction at https://github.com/brentp/smoove
# say that this command can be parallalized
# but smoove merge doesn't appear to have 
# a way to specify the number of threads
# I gave it 10 cpus to see if it will make use of them


# this will create ./merged.sites.vcf.gz
smoove merge --name merged -f $fasta --outdir $out_dir ${in_dir}/*.genotyped.vcf.gz





echo `date` finished



