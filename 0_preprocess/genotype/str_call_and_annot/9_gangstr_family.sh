#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=gangstr
#SBATCH -o logs/gangstr/submit_%a.log
#SBATCH -e logs/gangstr/submit_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=40G
#SBATCH --chdir /data5/16p12_WGS/structural_variants/gangstr
#SBATCH --nodelist=qingyu
#SBATCH --array=1


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

echo `date` $family
echo `date` $samples
echo `date` $bamfiles

/data5/software/GangSTR-2.5/build/GangSTR \
		--bam $bamfiles \
		--ref $fasta \
		--regions output/bed_non_ref_calls/$family.bed \
		--out $out_dir/$family \
		--include-ggl 

mv ${family}.vcf ${family}.bk.vcf
bgzip ${out_dir}/${family}.vcf
tabix -p vcf ${out_dir}/${family}.vcf.gz




echo `date` "done"

