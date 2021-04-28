#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=smoove
#SBATCH -o logs/smoove/submit_%a.log
#SBATCH -e logs/smoove/submit_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=5G
#SBATCH --workdir /data5/16p12_WGS/structural_variants/vcf_callers
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user awt5304@psu.edu
#SBATCH --nodelist=qingyu
#SBATCH --array=2-345%3

echo `date` starting job on $HOSTNAME as $USER

fasta=/data/bx_references/hg19/ucsc.hg19.fasta
out_dir=/data5/16p12_WGS/structural_variants/vcf_callers/output/smoove
tmp_dir=/data5/16p12_WGS/structural_variants/vcf_callers/tmp/smoove

bamfile=`cat bam.list | awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id'`
sample=`echo $bamfile | cut -d"/" -f5`

export PATH=$PATH:/data3/software/bin
export PATH=$PATH:/data3/software/lumpy-sv/bin

echo `date` $bamfile
echo `date` $sample

exclude_chroms=chr6_ssto_hap7,chr6_mcf_hap5,chr6_cox_hap2,chr6_mann_hap4,chr6_apd_hap1,chr6_qbl_hap6,chr6_dbb_hap3,chr17_ctg5_hap1,chr4_ctg9_hap1,chr1_gl000192_random,chrUn_gl000225,chr4_gl000194_random,chr4_gl000193_random,chr9_gl000200_random,chrUn_gl000222,chrUn_gl000212,chr7_gl000195_random,chrUn_gl000223,chrUn_gl000224,chrUn_gl000219,chr17_gl000205_random,chrUn_gl000215,chrUn_gl000216,chrUn_gl000217,chr9_gl000199_random,chrUn_gl000211,chrUn_gl000213,chrUn_gl000220,chrUn_gl000218,chr19_gl000209_random,chrUn_gl000221,chrUn_gl000214,chrUn_gl000228,chrUn_gl000227,chr1_gl000191_random,chr19_gl000208_random,chr9_gl000198_random,chr17_gl000204_random,chrUn_gl000233,chrUn_gl000237,chrUn_gl000230,chrUn_gl000242,chrUn_gl000243,chrUn_gl000241,chrUn_gl000236,chrUn_gl000240,chr17_gl000206_random,chrUn_gl000232,chrUn_gl000234,chr11_gl000202_random,chrUn_gl000238,chrUn_gl000244,chrUn_gl000248,chr8_gl000196_random,chrUn_gl000249,chrUn_gl000246,chr17_gl000203_random,chr8_gl000197_random,chrUn_gl000245,chrUn_gl000247,chr9_gl000201_random,chrUn_gl000235,chrUn_gl000239,chr21_gl000210_random,chrUn_gl000231,chrUn_gl000229,chrUn_gl000226,chr18_gl000207_random

smoove call --name $sample \
 --fasta $fasta \
 --outdir $tmp_dir \
 --excludechroms $exclude_chroms \
 $bamfile

echo `date` "done"

