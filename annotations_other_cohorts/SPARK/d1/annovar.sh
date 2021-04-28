#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=annovar
#SBATCH -o logs/annovar/submit_%a.log
#SBATCH -e logs/annovar/submit_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=2G
#SBATCH --chdir /data5/SPARK_WES/annotations
#SBATCH --array 0-9

# Note: don't need to do split alleles or left norm because the vcf files are already split and left normed

echo `date` started on $HOSTNAME

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

	vcf=vcfs/ref_filtered/${sample}.vcf.gz

	perl /data4/software/annovar/table_annovar.pl $vcf /data4/software/annovar/humandb/ \
	 -buildver hg38 \
	 -out vcfs/annovar/${sample} \
	 -remove \
	 -protocol refGene,knownGene \
	 -operation gx,gx \
	 -nastring . \
	 -thread 6 \
	 -vcfinput \
	 -xref /data4/software/annovar/humandb/gene_annotations.txt \
	 -arg '-hgvs','-hgvs'

	bgzip vcfs/annovar/$sample.hg38_multianno.vcf
	tabix -p vcf vcfs/annovar/$sample.hg38_multianno.vcf.gz


	# cleanup
	rm vcfs/annovar/$sample.avinput


done



echo `date` fin


