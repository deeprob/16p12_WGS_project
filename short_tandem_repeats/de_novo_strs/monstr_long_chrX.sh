#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=monstr_long
#SBATCH -o logs/monstr/submit_long_chrX_%a.log
#SBATCH -e logs/monstr/submit_long_chrX_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --chdir /data5/16p12_WGS/structural_variants/de_novo_strs
#SBATCH --nodelist=durga
#SBATCH --array=1-115

# https://www.nature.com/articles/s41586-020-03078-7#Sec6
# To additionally identify candidate expansions, we called MonSTR again on each family using the non-default parameter --naive-expansions-frr 3,8 



echo `date` starting job on $HOSTNAME as $USER

export PATH=$PATH:/data5/software/STRDenovoTools-1.0.0/build_durga

child=`cat children_of_trios.txt | awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id'`

ped=peds/$child.trio.ped
vcf=vcfs/chrX/$child.chrX.vcf.gz


out=output/monstr/long/${child}_chrX

echo `date` $vcf

MonSTR \
    --strvcf $vcf \
    --fam $ped \
    --gangstr \
    --out $out \
    --output-all-loci \
    --max-num-alleles 100 \
    --naive-expansions-frr 3,8 \
    --include-invariant \
	--require-all-children \
	--naive \
    --chrX





echo `date` "done"


