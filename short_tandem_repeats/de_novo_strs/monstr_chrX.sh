#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=monstr
#SBATCH -o logs/monstr/submit_chrX_%a.log
#SBATCH -e logs/monstr/submit_chrX_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --chdir /data5/16p12_WGS/structural_variants/de_novo_strs
#SBATCH --nodelist=durga
#SBATCH --array 34,21,37,35,49

# https://www.nature.com/articles/s41586-020-03078-7#Sec6
# MonSTR v1.0.0 was called separately on each family after applying call-level and locus-level genotype filters described above. MonSTR was called with non-default parameters --max-num-alleles 100 --include-invariant --gangstr --require-all-children --output-all-loci --min-num-encl-child 3 --max-perc-encl-parent 0.05 --min-encl-match 0.9 --min-total-encl 10 --posterior-threshold 0.5. Autosomes were run with the --default-prior -3 and chromosome X was run with the --naive option. 


echo `date` starting job on $HOSTNAME as $USER

export PATH=$PATH:/data5/software/STRDenovoTools-1.0.0/build_durga

child=`cat children_of_trios.txt | awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id'`

ped=peds/$child.trio.ped
vcf=vcfs/chrX/$child.chrX.vcf.gz

out=output/monstr/${child}.chrX

echo `date` $vcf

MonSTR \
    --strvcf $vcf \
    --fam $ped \
    --gangstr \
    --out $out \
	--max-num-alleles 100 \
	--include-invariant \
	--require-all-children \
	--output-all-loci \
	--min-num-encl-child 3 \
	--max-perc-encl-parent 0.05 \
	--min-encl-match 0.9 \
	--min-total-encl 10 \
	--posterior-threshold 0.5 \
	--naive \
	--chrX



echo `date` "done"





