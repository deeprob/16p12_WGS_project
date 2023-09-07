#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=filter1
#SBATCH -o logs/2_filter1/batch_%a.log
#SBATCH -e logs/2_filter1/batch_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=2G
#SBATCH --chdir /data5/16p12_WGS/annotations/2022_01_27/
#SBATCH --array 2-346%100
#SBATCH --nodelist qingyu

# need to filter for alt
# norm split
# and norm left
# QUAL>=50
# depth>=8
# alternative depth>0
# alternative allele fraction >=0.25 and <=0.75
# and QUAL/(alternative depth)>=1.5


# Note: maximum slurm array job size is 1000 

echo `date` starting job on $HOSTNAME


FASTA=/data/bx_references/hg19/ucsc.hg19.fasta

SAMPLE=`cat samples.csv | head -n $SLURM_ARRAY_TASK_ID | tail -n1 | cut -f2 -d,`
INPUT_VCF=`cat samples.csv | head -n $SLURM_ARRAY_TASK_ID | tail -n1 | cut -f1 -d,`
OUTPUT_FILE=/data5/16p12_WGS/annotations/2022_01_27/vcfs/filter1/$SAMPLE.vcf.gz

echo $SAMPLE
echo input file: $INPUT_VCF
echo output file: $OUTPUT_FILE

export PATH=/data5/anastasia/sw/bcftools-1.12:$PATH


# bcftools norm -m-both $INPUT_VCF | bcftools view -i 'GT="alt" & QUAL>=50 & FORMAT/DP>=8 & FORMAT/AD[:1]>0 & (FORMAT/AD[:1])/(FORMAT/DP)>=0.25 & (FORMAT/AD[:1])/(FORMAT/DP)<=0.75 & QUAL/(FORMAT/AD[:1])>=1.5' |  bcftools norm -f ${FASTA} | bgzip > ${OUTPUT_FILE}
bcftools norm -m-both $INPUT_VCF | bcftools view -i 'GT="alt" & QUAL>=50 & FORMAT/DP>=8 & FORMAT/AD[:1]>0 & (FORMAT/AD[:1])/(FORMAT/DP)>=0.25 & QUAL/(FORMAT/AD[:1])>=1.5' | bcftools view -i '(FORMAT/AD[:1])/(FORMAT/DP)<=0.75 | (FORMAT/AD[:1])/(FORMAT/DP)>=0.9' |  bcftools norm -f ${FASTA} | bgzip > ${OUTPUT_FILE}


echo `date` finished












