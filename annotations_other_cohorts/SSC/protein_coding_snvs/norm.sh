#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=norm
#SBATCH -o logs/norm_%a.log
#SBATCH -e logs/norm_%a.err
#SBATCH --cpus-per-task=2
#SBATCH --time=240:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --chdir /data5/SSC_WGS/SFARI_SSC_WGS_2/annotations
#SBATCH --array 1-20,23
#SBATCH --nodelist durga


echo `date` started on $HOSTNAME


export PATH=/data5/anastasia/sw/bcftools-1.12/build/bin/:$PATH

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


vcfdir=/data5/SSC_WGS/SFARI_SSC_WGS_2/Project_SSC_9209samples.JGvariants.2019-06-21/

vcf=${vcfdir}/CCDG_9000JG_B01_GRM_WGS_2019-03-21_chr${a}.recalibrated_variants.flagged.vcf.gz
out=vcfs/norm/chr${a}.vcf.gz

fasta=../documentation/GRCh38_full_analysis_set_plus_decoy_hla.fa


echo `date` working on chr${a}

bcftools norm -f $fasta -m-both $vcf | bgzip --threads 4 > $out
tabix -p vcf $out

echo `date` fin


