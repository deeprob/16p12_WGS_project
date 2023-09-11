#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=get_calls
#SBATCH -o logs/get_calls/submit_%a.log
#SBATCH -e logs/get_calls/submit_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=2G
#SBATCH --chdir /data5/16p12_WGS/structural_variants/gangstr
#SBATCH --array 1-129%10

echo `date` starting job on $HOSTNAME as $USER

family=`cat families_list.txt | awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id'`


vcf=output/vcfs_filtered/${family}.vcf.gz
bed=../vcf_callers/hg19_ver13_1.bed
out=output/vcfs_noref/${family}.vcf.gz

bcftools view -i 'GT[*]="alt"' $vcf | bgzip > $out


tabix -p vcf $out






echo `date` "done"
