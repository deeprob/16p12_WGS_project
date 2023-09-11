#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=dumpstr
#SBATCH -o logs/dumpstr/submit_%a.log
#SBATCH -e logs/dumpstr/submit_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=2G
#SBATCH --workdir /data5/16p12_WGS/structural_variants/vcf_callers
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user awt5304@psu.edu
#SBATCH --nodelist=qingyu
#SBATCH --array=1-345

echo `date` starting job on $HOSTNAME as $USER

export SINGULARITY_CACHEDIR=/data5/16p12_WGS/structural_variants/vcf_callers/cache


bamfile=`cat bam.list | awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id'`
sample=`echo $bamfile | cut -d"/" -f5`
vcf=output/gangstr/gangstr.${sample}.vcf.gz


echo `date` $vcf
echo `date` $sample

singularity exec \
	-B /data5 \
	/data5/software/str-toolkit_latest.sif \
	dumpSTR \
	--verbose \
	--vcf $vcf \
	--out tmp/dumpstr/$sample \
	--min-call-DP 20 \
	--max-call-DP 1000 \
	--filter-spanbound-only \
	--filter-badCI \
	--min-call-Q 0.9

cp tmp/dumpstr/${sample}.vcf output/dumpstr/dumpstr.${sample}.vcf
bgzip output/dumpstr/dumpstr.${sample}.vcf
tabix -p vcf output/dumpstr/dumpstr.${sample}.vcf.gz


echo `date` "done"




