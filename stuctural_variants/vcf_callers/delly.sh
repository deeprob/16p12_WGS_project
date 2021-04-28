#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=delly
#SBATCH -o logs/delly/submit_%a.log
#SBATCH -e logs/delly/submit_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=10G
#SBATCH --workdir /data5/16p12_WGS/structural_variants/vcf_callers
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user awt5304@psu.edu
#SBATCH --nodelist=qingyu
#SBATCH --array=292-345

echo `date` starting job on $HOSTNAME as $USER

fasta=/data/bx_references/hg19/ucsc.hg19.fasta
out_dir=/data5/16p12_WGS/structural_variants/vcf_callers/output/delly
tmp_dir=/data5/16p12_WGS/structural_variants/vcf_callers/tmp/delly

bamfile=`cat bam.list | awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id'`
sample=`echo $bamfile | cut -d"/" -f5`

echo `date` $bamfile
echo `date` $sample

delly=/data3/software/delly_v0.8.2_linux_x86_64bit

$delly call \
	-o ${tmp_dir}/${sample}.bcf \
	-g $fasta \
	-x exclude.bed \
	$bamfile

bcftools view ${tmp_dir}/${sample}.bcf > ${out_dir}/${sample}.delly.vcf


echo `date` "done"


