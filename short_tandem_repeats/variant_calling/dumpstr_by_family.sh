#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=dumpstr
#SBATCH -o logs/dumpstr_by_family/submit_%a.log
#SBATCH -e logs/dumpstr_by_family/submit_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=2G
#SBATCH --chdir /data5/16p12_WGS/structural_variants/vcf_callers
#SBATCH --nodelist=durga
#SBATCH --array=1-129

echo `date` starting job on $HOSTNAME as $USER

export SINGULARITY_CACHEDIR=/data5/16p12_WGS/structural_variants/vcf_callers/durga_cache

out_dir=/data5/16p12_WGS/structural_variants/vcf_callers/output/dumpstr_by_family
tmp_dir=/data5/16p12_WGS/structural_variants/vcf_callers/tmp/dumpstr_by_family

family=`cat families_list.txt | awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id'`
vcf=/data5/16p12_WGS/structural_variants/vcf_callers/output/gangstr_by_family/gangstr.${family}.vcf.gz

echo `date` $family
echo `date` $vcf

singularity exec \
	-B /data5 \
	durga_cache/str-toolkit.simg \
	dumpSTR \
	--verbose \
	--vcf $vcf \
	--out ${tmp_dir}/$family \
	--min-call-DP 20 \
	--max-call-DP 1000 \
	--filter-spanbound-only \
	--filter-badCI \
	--readlen 150
	
# Note: these are the same parameters as Gymrek's paper (excluding --require-support 2 which results in empty VCF)
# see https://www.nature.com/articles/s41586-020-03078-7#Sec6

cp ${tmp_dir}/${family}.vcf ${out_dir}/dumpstr.${family}.vcf
bgzip ${out_dir}/dumpstr.${family}.vcf
tabix -p vcf ${out_dir}/dumpstr.${family}.vcf.gz

echo `date` "done"




