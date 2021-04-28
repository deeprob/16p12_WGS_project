#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=ref_fil
#SBATCH -o logs/ref_filtered_%a.log
#SBATCH -e logs/ref_filtered_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=2G
#SBATCH --chdir /data5/SPARK_WES/annotations/release_june2020
#SBATCH --array 0-9


# Removes vcf records where both alleles are REF to significantly reduce file size

# Note: maximum slurm array job size is 1000 and the number of samples is 15995

echo `date` starting job on $HOSTNAME as $USER

# annotations for samples ending with $a
a=$SLURM_ARRAY_TASK_ID
echo "--annotations for samples ending with $a"

# get samples to annotate
all_samples=`cat samples_june2020.txt | grep "${a}$"`

# num samples to annotate
num_samples=`cat samples_june2020.txt | grep "${a}$" | wc -l`
echo "--number of annotations: $num_samples"


# look through samples
for i in `seq 1 $num_samples`
do

sample=`echo $all_samples | cut -f $i -d' '`
echo "--working on ${i}/${num_samples} $sample"

# get last digit
last_digit=`echo "${sample: -1}"`

# vcf file
vcf=/data5/SPARK_WES/Variants/GATK/gvcf_individual/$last_digit/${sample}.gvcf.gz

# outfile
out=vcfs/ref_filtered/${sample}.vcf.gz

if [ -e $out ]
then
	# do nothing
	echo 'already done'
else
	echo `date` $vcf
	echo `date` $out

	python3 ../filter_ref.py $vcf | bgzip > $out
fi


done

echo `date` finished

