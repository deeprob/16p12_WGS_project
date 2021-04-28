#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=gangstr_by_family
#SBATCH -o logs/gangstr_by_family_chrX/submit_%a.log
#SBATCH -e logs/gangstr_by_family_chrX/submit_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=40G
#SBATCH --chdir /data5/16p12_WGS/structural_variants/vcf_callers
#SBATCH --nodelist=qingyu
#SBATCH --array=1-129

# https://www.nature.com/articles/s41586-020-03078-7#Sec6
# From Gymrek autism paper, "Chromosome X TRs were genotyped using GangSTR v2.4.4 with additional options --bam-samps and --samp-sex to interpret sample sex for chromosome X. A separate GangSTR job was run for each family on each chromosome resulting in separate VCF files for each.""

# Here, we use the old option --ploidy rather than the new option --samp-sex because samp-sex doesn't work well with dumpSTR
# --samp-sex outputs only 1 GT field per male chromosome X, which dumpSTR doesn't like

echo `date` starting job on $HOSTNAME as $USER

fasta=/data/bx_references/hg19/ucsc.hg19.fasta
out_dir=/data5/16p12_WGS/structural_variants/vcf_callers/output/gangstr_by_family_chrX
tmp_dir=/data5/16p12_WGS/structural_variants/vcf_callers/tmp/gangstr_by_family_chrX

family=`cat families_list.txt | awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id'`
samples=`cat all_batches.ped | grep $family | cut -f2`

bamfiles=
for samp in $samples
do
append_bam=`cat bam.list | grep $samp`
bamfiles=`echo $bamfiles,$append_bam`
done
bamfiles=`echo ${bamfiles:1}`

samplesex=
for samp in $samples
do
append_sex=`cat sample_sex.txt | grep $samp | cut -f2`
samplesex=`echo $samplesex,$append_sex`
done
samplesex=`echo ${samplesex:1}`

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
