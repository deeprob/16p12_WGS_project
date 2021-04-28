#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=mergestr
#SBATCH -o logs/mergestr_by_family/submit.log
#SBATCH -e logs/mergestr_by_family/submit.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=2G
#SBATCH --chdir /data5/16p12_WGS/structural_variants/vcf_callers
#SBATCH --nodelist=durga

echo `date` starting job on $HOSTNAME as $USER

export SINGULARITY_CACHEDIR=/data5/16p12_WGS/structural_variants/vcf_callers/durga_cache

out_dir=/data5/16p12_WGS/structural_variants/vcf_callers/output/mergestr_by_family

families=`cat families_list.txt`

# get all vcfs ============
vcffiles=
for fam in $families
do
append_vcf=output/dumpstr_by_family/dumpstr.${fam}.vcf.gz
vcffiles=`echo $vcffiles,$append_vcf`
done
vcffiles=`echo ${vcffiles:1}`
# =========================

echo `date` $families
echo `date` $vcffiles

singularity exec \
	-B /data5 \
	durga_cache/str-toolkit.simg \
	mergeSTR \
	--vcfs $vcffiles \
	--out $out_dir/mergestr

echo `date` "done"




