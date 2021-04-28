#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=get_calls
#SBATCH -o logs/bed_non_ref_calls/submit_%a.log
#SBATCH -e logs/bed_non_ref_calls/submit_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=2G
#SBATCH --chdir /data5/16p12_WGS/structural_variants/gangstr
#SBATCH --array 1-129%10

echo `date` starting job on $HOSTNAME as $USER

family=`cat families_list.txt | awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id'`


vcf=output/vcfs_noref/${family}.vcf.gz
out=output/bed_non_ref_calls/${family}.bed


vcf-query $vcf -f '%CHROM\t%POS\t%INFO/END\t%INFO/PERIOD\t%INFO/RU\n'  > $out

# vcf-query $vcf -f '%CHROM\t%POS\t%INFO/END\t%INFO/PERIOD\t%INFO/RU\n' | awk '$5 = toupper($5)' > $out



echo `date` "done"
