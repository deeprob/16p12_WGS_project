#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=gangstr
#SBATCH -o logs/gangstr/submit_PSU008.log
#SBATCH -e logs/gangstr/submit_PSU008.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=40G
#SBATCH --chdir /data5/16p12_WGS/structural_variants/gangstr
#SBATCH --nodelist=qingyu


echo `date` starting job on $HOSTNAME as $USER


fasta=/data/bx_references/hg19/ucsc.hg19.fasta
out_dir=output/gangstr

family=PSU008
samples="SG047 SG050"

bamfiles="/data2/16p12_WGS/batch2_0517/SG047/SG047_rmdup_BQSR.bam,/data3/16p12_WGS/batch3_0618/SG050/SG050_rmdup_BQSR.bam"

echo `date` $family
echo `date` $samples
echo `date` $bamfiles

/data5/software/GangSTR-2.5/build/GangSTR \
		--bam $bamfiles \
		--ref $fasta \
		--regions output/bed_non_ref_calls/$family.bed \
		--out $out_dir/$family \
		--include-ggl 

bgzip ${out_dir}/${family}.vcf
tabix -p vcf ${out_dir}/${family}.vcf.gz




echo `date` "done"

