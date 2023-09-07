#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=gnomad
#SBATCH -o logs/6_filter_pop_freq/batch_%a.log
#SBATCH -e logs/6_filter_pop_freq/batch_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=10G
#SBATCH --chdir /data5/16p12_WGS/annotations/2022_02_04/
#SBATCH --array 2-346%100
#SBATCH --nodelist qingyu

echo `date` starting job on $HOSTNAME


export PATH=/data5/anastasia/sw/bcftools-1.12:$PATH


SAMPLE=`cat samples.csv | head -n $SLURM_ARRAY_TASK_ID | tail -n1 | cut -f2 -d,`
INPUT_VCF=/data5/16p12_WGS/annotations/2022_02_04/exonic_vcfs/gnomad/$SAMPLE.vcf.gz
OUTPUT_VCF=/data5/16p12_WGS/annotations/2022_02_04/exonic_vcfs/filter2/$SAMPLE.vcf.gz


echo $SAMPLE
echo input file: $INPUT_VCF
echo output file: $OUTPUT_VCF



bcftools view -i 'gnomad_exome_AF<=0.001 | gnomad_exome_AF="."' $INPUT_VCF | bcftools view -i 'gnomad_genome_AF<=0.001 | gnomad_genome_AF="."' | bgzip > $OUTPUT_VCF





 
echo `date` finished

 
 
 