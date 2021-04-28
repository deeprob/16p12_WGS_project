#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=gangstr_chrX
#SBATCH -o logs/gangstr/submit_chrX_%a.log
#SBATCH -e logs/gangstr/submit_chrX_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=40G
#SBATCH --chdir /data5/16p12_WGS/structural_variants/gangstr
#SBATCH --nodelist=qingyu
#SBATCH --array=11-129%10


echo `date` starting job on $HOSTNAME as $USER


fasta=/data/bx_references/hg19/ucsc.hg19.fasta
out_dir=output/gangstr

family=`cat families_list.txt | awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id'`
samples=`cat all_batches.ped | grep $family | cut -f2`

bamfiles=
for samp in $samples
do
append_bam=`cat ../vcf_callers/bam.list | grep $samp`
bamfiles=`echo $bamfiles,$append_bam`
done
bamfiles=`echo ${bamfiles:1}`

samplesex=
for samp in $samples
do
append_sex=`cat ../vcf_callers/sample_sex.txt | grep $samp | cut -f3`
samplesex=`echo $samplesex,$append_sex`
done
samplesex=`echo ${samplesex:1}`

bamsamps=`echo $samples | sed 's/\s\+/,/g'`

echo `date` $family
echo `date` $samples
echo `date` $bamfiles
echo `date` $samplesex
echo `date` $bamsamps

/data5/software/GangSTR-2.5/build/GangSTR \
		--bam $bamfiles \
 		--ref $fasta \
 		--regions output/bed_non_ref_calls/$family.chrX.bed  \
 		--out $out_dir/$family.chrX \
 		--include-ggl \
 		--chrom chrX \
 		--samp-sex $samplesex \
		--bam-samps $bamsamps

bgzip ${out_dir}/${family}.chrX.vcf
tabix -p vcf ${out_dir}/${family}.chrX.vcf.gz




echo `date` "done"


