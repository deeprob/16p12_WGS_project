#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=gangstr_by_family
#SBATCH -o logs/gangstr_by_family_chrX/submit_PSU008.log
#SBATCH -e logs/gangstr_by_family_chrX/submit_PSU008.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=40G
#SBATCH --chdir /data5/16p12_WGS/structural_variants/vcf_callers
#SBATCH --nodelist=qingyu

# Same as gangstr_by_family_chromX.sh
# except this is only for PSU008

echo `date` starting job on $HOSTNAME as $USER

fasta=/data/bx_references/hg19/ucsc.hg19.fasta
out_dir=/data5/16p12_WGS/structural_variants/vcf_callers/output/gangstr_by_family_chrX
tmp_dir=/data5/16p12_WGS/structural_variants/vcf_callers/tmp/gangstr_by_family_chrX

family=PSU008
samples="SG047 SG050"

bamfiles=
for samp in $samples
do
append_bam=`cat bam.list | grep $samp`
bamfiles=`echo $bamfiles,$append_bam`
done
bamfiles=`echo ${bamfiles:1}`

samplesex="2,2"

echo `date` $family
echo `date` $samples
echo `date` $bamfiles
echo `date` $samplesex

singularity exec -B /data/bx_references/hg19,/data3/16p12_WGS,/data4/16p12_WGS,/data2/16p12_WGS,/data5 \
	durga_cache/str-toolkit.simg  \
	GangSTR \
		--bam $bamfiles \
		--ref $fasta \
		--regions hg19_ver13_1.bed \
		--out $tmp_dir/$family \
		--include-ggl \
		--chrom chrX \
		--ploidy $samplesex

		

cp ${tmp_dir}/${family}.vcf ${out_dir}/gangstr.${family}.vcf
bgzip ${out_dir}/gangstr.${family}.vcf
tabix -p vcf ${out_dir}/gangstr.${family}.vcf.gz




echo `date` "done"
