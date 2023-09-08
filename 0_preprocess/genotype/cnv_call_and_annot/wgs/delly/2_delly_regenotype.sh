#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=delly
#SBATCH -o logs/2_regenotype/%a.log
#SBATCH -e logs/2_regenotype/%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=5G
#SBATCH --chdir /data5/16p12_WGS/structural_variants/vcf_callers_2022_02_07/delly
#SBATCH --array=1-345%100



echo `date` starting job on $HOSTNAME as $USER


bamfile=`head -n $SLURM_ARRAY_TASK_ID /data5/16p12_WGS/structural_variants/vcf_callers/bam.list | tail -n 1`
sample=`echo $bamfile | cut -d"/" -f5`

delly=/data3/software/delly_v0.8.2_linux_x86_64bit
fasta=/data5/bx_reference/hg19/ucsc.hg19.fasta


# delly call -g hg19.fa -v sites.bcf -o s1.geno.bcf -x hg19.excl s1.bam
$delly call -g $fasta -v /data5/16p12_WGS/structural_variants/vcf_callers_2022_02_07/delly/vcfs/1_merge/merged.bcf -o /data5/16p12_WGS/structural_variants/vcf_callers_2022_02_07/delly/vcfs/2_regenotype/$sample.bcf $bamfile


echo `date` finished




