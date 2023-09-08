#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=smoove
#SBATCH -o logs/3_regenotype/%a.log
#SBATCH -o logs/3_regenotype/%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=5G
#SBATCH --chdir /data5/16p12_WGS/structural_variants/vcf_callers_2022_02_07/smoove_lumpy
#SBATCH --array=1-345%30
#SBATCH --nodelist qingyu



echo `date` starting job on $HOSTNAME as $USER


fasta=/data5/bx_reference/hg19/ucsc.hg19.fasta
in_file=/data5/16p12_WGS/structural_variants/vcf_callers_2022_02_07/smoove_lumpy/vcfs/2_union/merged.sites.vcf.gz
out_dir=/data5/16p12_WGS/structural_variants/vcf_callers_2022_02_07/smoove_lumpy/vcfs/3_regenotype/


bamfile=`head -n $SLURM_ARRAY_TASK_ID /data5/16p12_WGS/structural_variants/vcf_callers/bam.list | tail -n 1`
sample=`echo $bamfile | cut -d"/" -f5`


export TMPDIR=/data5/16p12_WGS/structural_variants/vcf_callers_2022_02_07/smoove_lumpy/tmp
export PATH=/data5/anastasia/sw/bcftools-1.12:$PATH
export PATH=/data3/software/bin:$PATH
export PATH=/data5/software/lumpy-sv/bin/:$PATH
export PATH=$PATH:/data5/software/


# Instructions at https://github.com/brentp/smoove
# run on each sample at a time
# can remerge into a single cohort vcf later if needed


smoove genotype -d -x -p 20 --name $sample --outdir $out_dir --fasta $fasta --vcf $in_file $bamfile





echo `date` "done"





