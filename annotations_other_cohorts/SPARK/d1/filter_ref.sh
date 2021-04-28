#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=ref_fil
#SBATCH -o logs/ref_filtered/submit_rest.log
#SBATCH -e logs/ref_filtered/submit_rest.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=2G
#SBATCH --chdir /data5/SPARK_WES/annotations

# Removes vcf records where both alleles are REF to significantly reduce file size

# Note: maximum slurm array job size is 1000 and the number of samples is 27290
# This script is quick so for samples 1001-27290 I'll use a loop

echo `date` starting job on $HOSTNAME as $USER


for i in `seq 1001 27290`
do

sample=`cat samples_d1.txt | awk -v task_id=$i 'NR==task_id'`

# get last digit
last_digit=`echo "${sample: -1}"`

# vcf file
vcf=/data4/SPARK/Data1_Regeneron/wes.gVCF/$last_digit/${sample}.acmg56excluded.g.vcf.gz

# outfile
out=vcfs/ref_filtered/${sample}.vcf.gz

echo `date` $i $sample $last_digit
echo `date` $vcf
echo `date` $out

python3 filter_ref.py $vcf | bgzip > $out

done

echo `date` finished
