#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=gnomad
#SBATCH -o logs/gnomad_%a.log
#SBATCH -e logs/gnomad_%a.err
#SBATCH --cpus-per-task=2
#SBATCH --time=240:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --chdir /data5/SSC_WGS/SFARI_SSC_WGS_2/annotations
#SBATCH --array 1-24

echo `date` started on $HOSTNAME

export PATH=$PATH:/data5/software

vcfdir=/data5/SSC_WGS/SFARI_SSC_WGS_2/Project_SSC_9209samples.JGvariants.2019-06-21/

# annotations for chroms ending with $a
a=$SLURM_ARRAY_TASK_ID

# if s > 22 set x,y,m chromosomes
if [[ $a -eq 23 ]]
then
  a=X
elif [[ $a -eq 24 ]]
then
  a=Y
fi

vcf=${vcfdir}/CCDG_9000JG_B01_GRM_WGS_2019-03-21_chr${a}.recalibrated_variants.annotated.vcf.gz
out=vcfs/gnomad/chr${a}.vcf.gz

echo `date` working on chr${a}
echo `date` invcf=$vcf
echo `date` outvcf=$out

slivar expr --vcf $vcf \
	-g /data5/anastasia/slivar/hg38/annotations/gnomad.zip \
	| bgzip \
	> $out

tabix -p vcf $out









echo `date` fin




