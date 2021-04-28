#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=MPC
#SBATCH -o logs/MPC/submit_%a.log
#SBATCH -e logs/MPC/submit_%a.err
#SBATCH --cpus-per-task=2
#SBATCH --time=240:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --chdir /data5/SPARK_WES/annotations
#SBATCH --array 0-9


echo `date` started on $HOSTNAME

export PATH=$PATH:/data5/software

# annotations for samples ending with $a
a=$SLURM_ARRAY_TASK_ID
echo "--annotations for samples ending with $a"

# get samples to annotate
all_samples=`cat samples_d1.txt | grep "${a}$"`

# num samples to annotate
num_samples=`cat samples_d1.txt | grep "${a}$" | wc -l`
echo "--number of annotations: $num_samples"

# look through samples
for i in `seq 1 $num_samples`
do

	sample=`echo $all_samples | cut -f $i -d' '`
	echo "--working on ${i}/${num_samples} $sample"

	vcf=vcfs/CADD/${sample}.vcf.gz
	out=vcfs/MPC/${sample}.vcf.gz

	if [ -e $out.tbi ]
	then
		# do nothing
		echo 'already done'
	else
		slivar expr --vcf $vcf \
	        -g /data5/anastasia/liftover/mpc_hg38_gnotate.zip \
    		| bgzip \
    		> $out

		tabix -p vcf $out
	fi



done






echo `date` fin









