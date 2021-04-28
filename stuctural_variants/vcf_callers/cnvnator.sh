#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=cnvnator
#SBATCH -o logs/cnvnator/submit_%a.log
#SBATCH -e logs/cnvnator/submit_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=3G
#SBATCH --workdir /data5/16p12_WGS/structural_variants/vcf_callers
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user awt5304@psu.edu
#SBATCH --nodelist=qingyu
#SBATCH --array=1-243%5

echo `date` starting job on $HOSTNAME as $USER

source /data5/software/root/bin/thisroot.sh
cnvnator_dir=/data5/software/CNVnator_v0.4.1/src/
export PATH=$PATH:$cnvnator_dir

out_dir=/data5/16p12_WGS/structural_variants/vcf_callers/output/cnvnator
tmp_dir=/data5/16p12_WGS/structural_variants/vcf_callers/tmp/cnvnator

bamfile=`cat bam.list | awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id'`
sample=`echo $bamfile | cut -d"/" -f5`

echo `date` $bamfile
echo `date` $sample

cnvnator -root ${tmp_dir}/${sample}.root -tree $bamfile -chrom $(seq -f 'chr%g' 1 22) chrX chrY chrM # 14 gb
cnvnator -root ${tmp_dir}/${sample}.root -his 200 -chrom $(seq -f 'chr%g' 1 22) chrX chrY chrM -d $cnvnator_dir/../fastas # 2 gb
cnvnator -root ${tmp_dir}/${sample}.root -stat 200 # 4gb
cnvnator -root ${tmp_dir}/${sample}.root -partition 200 # 12214% CPU # 2 gb
cnvnator -root ${tmp_dir}/${sample}.root -call 200 > output/cnvnator/bin200/${sample}.cnvnator.vcf # 4 gb

echo `date` "done"



