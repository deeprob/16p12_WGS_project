#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=filter
#SBATCH -o logs/filter_vcf/submit_%a.log
#SBATCH -e logs/filter_vcf/submit_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=2G
#SBATCH --chdir /data5/16p12_WGS/structural_variants/gangstr
#SBATCH --array 1-129%10

echo `date` starting job on $HOSTNAME as $USER

family=`cat families_list.txt | awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id'`

vcf=/data5/16p12_WGS/structural_variants/vcf_callers/output/dumpstr_by_family/dumpstr.${family}.vcf.gz
bed=output/merged_loci_remove.bed
out_prefix=output/vcfs_filtered/${family}

vcftools --exclude-bed $bed --gzvcf $vcf --recode --recode-INFO-all --out $out_prefix

cat $out_prefix.recode.vcf | vcf-sort | bgzip -c > $out_prefix.vcf.gz

tabix -p vcf $out_prefix.vcf.gz



echo `date` "done"




