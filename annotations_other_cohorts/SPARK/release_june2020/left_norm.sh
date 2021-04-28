#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=left_norm
#SBATCH -o logs/left_norm_%a.log
#SBATCH -e logs/left_norm_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=2G
#SBATCH --chdir /data5/SPARK_WES/annotations/release_june2020
#SBATCH --array 0-9


echo `date` starting job on $HOSTNAME as $USER


export PATH=/data5/anastasia/sw/bcftools-1.12/build/bin/:$PATH

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

if [ -e vcfs/left_norm/${sample}.vcf.gz ]
then
# do nothing
echo 'already done'
else

# trim alt alleles
in=vcfs/ref_filtered/${sample}.vcf.gz
out=vcfs/alt_trimmed/${sample}.vcf.gz
bcftools view --trim-alt-alleles $in | bgzip > $out

# left norm
in=vcfs/alt_trimmed/${sample}.vcf.gz
fasta=/data5/SPARK_WES/Resources/genome.fa
out=vcfs/left_norm/${sample}.vcf.gz
bcftools norm -f $fasta -m-both $in | bgzip > $out

fi
done

echo `date` finished
