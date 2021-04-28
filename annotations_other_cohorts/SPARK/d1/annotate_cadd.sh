#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=CADD
#SBATCH -o logs/CADD/submit_%a.log
#SBATCH -e logs/CADD/submit_%a.err
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

# LARGE Function that annotates a sample
# skip below for the for loop
#=========================================
function annotate_cadd {

sample=$1
echo $sample

vcf=vcfs/gnomad/${sample}.vcf.gz

# annotate each chrom

for chrom in `seq 6 22`
do

	out=vcfs/by_chrom/${sample}.${chrom}.cadd.vcf.gz

	slivar expr --vcf $vcf \
	    		-g /data5/anastasia/slivar/hg38/annotations/by_chrom/cadd.chr${chrom}.zip \
	    		--region $chrom \
	    		| bgzip \
	    		> $out

done

chrom=1
endpoint=249250621
midpoint=124625310
midpointp1=124625311

out=vcfs/by_chrom/${sample}.${chrom}.part1.cadd.vcf.gz
slivar expr --vcf $vcf \
    		-g /data5/anastasia/slivar/hg38/annotations/by_chrom/cadd.chr${chrom}.part1.zip \
    		--region ${chrom}:0-${midpoint} \
    		| bgzip \
    		> $out

out=vcfs/by_chrom/${sample}.${chrom}.part2.cadd.vcf.gz
slivar expr --vcf $vcf \
    		-g /data5/anastasia/slivar/hg38/annotations/by_chrom/cadd.chr${chrom}.part1.zip \
    		--region ${chrom}:${midpointp1}-${endpoint} \
    		| bgzip \
    		> $out

chrom=2
endpoint=243199373
midpoint=121599686
midpointp1=121599687

out=vcfs/by_chrom/${sample}.${chrom}.part1.cadd.vcf.gz
slivar expr --vcf $vcf \
    		-g /data5/anastasia/slivar/hg38/annotations/by_chrom/cadd.chr${chrom}.part1.zip \
    		--region ${chrom}:0-${midpoint} \
    		| bgzip \
    		> $out

out=vcfs/by_chrom/${sample}.${chrom}.part2.cadd.vcf.gz
slivar expr --vcf $vcf \
    		-g /data5/anastasia/slivar/hg38/annotations/by_chrom/cadd.chr${chrom}.part1.zip \
    		--region ${chrom}:${midpointp1}-${endpoint} \
    		| bgzip \
    		> $out

chrom=3
endpoint=198295559
midpoint=99011215
midpointp1=99011216

out=vcfs/by_chrom/${sample}.${chrom}.part1.cadd.vcf.gz
slivar expr --vcf $vcf \
    		-g /data5/anastasia/slivar/hg38/annotations/by_chrom/cadd.chr${chrom}.part1.zip \
    		--region ${chrom}:0-${midpoint} \
    		| bgzip \
    		> $out

out=vcfs/by_chrom/${sample}.${chrom}.part2.cadd.vcf.gz
slivar expr --vcf $vcf \
    		-g /data5/anastasia/slivar/hg38/annotations/by_chrom/cadd.chr${chrom}.part1.zip \
    		--region ${chrom}:${midpointp1}-${endpoint} \
    		| bgzip \
    		> $out

chrom=4
endpoint=191154276
midpoint=95577138
midpointp1=95577139

out=vcfs/by_chrom/${sample}.${chrom}.part1.cadd.vcf.gz
slivar expr --vcf $vcf \
    		-g /data5/anastasia/slivar/hg38/annotations/by_chrom/cadd.chr${chrom}.part1.zip \
    		--region ${chrom}:0-${midpoint} \
    		| bgzip \
    		> $out

out=vcfs/by_chrom/${sample}.${chrom}.part2.cadd.vcf.gz
slivar expr --vcf $vcf \
    		-g /data5/anastasia/slivar/hg38/annotations/by_chrom/cadd.chr${chrom}.part1.zip \
    		--region ${chrom}:${midpointp1}-${endpoint} \
    		| bgzip \
    		> $out


chrom=5
endpoint=181538259
midpoint=90769128
midpointp1=90769129

out=vcfs/by_chrom/${sample}.${chrom}.part1.cadd.vcf.gz
slivar expr --vcf $vcf \
    		-g /data5/anastasia/slivar/hg38/annotations/by_chrom/cadd.chr${chrom}.part1.zip \
    		--region ${chrom}:0-${midpoint} \
    		| bgzip \
    		> $out

out=vcfs/by_chrom/${sample}.${chrom}.part2.cadd.vcf.gz
slivar expr --vcf $vcf \
    		-g /data5/anastasia/slivar/hg38/annotations/by_chrom/cadd.chr${chrom}.part1.zip \
    		--region ${chrom}:${midpointp1}-${endpoint} \
    		| bgzip \
    		> $out

chrom=X
out=vcfs/by_chrom/${sample}.${chrom}.cadd.vcf.gz
slivar expr --vcf $vcf \
    		-g /data5/anastasia/slivar/hg38/annotations/by_chrom/cadd.chr${chrom}.zip \
    		--region $chrom \
    		| bgzip \
    		> $out

chrom=Y
out=vcfs/by_chrom/${sample}.${chrom}.cadd.vcf.gz
slivar expr --vcf $vcf \
    		-g /data5/anastasia/slivar/hg38/annotations/by_chrom/cadd.chr${chrom}.zip \
    		--region $chrom \
    		| bgzip \
    		> $out

/data/software/bcftools/bcftools concat vcfs/by_chrom/${sample}.1.part1.cadd.vcf.gz \
vcfs/by_chrom/${sample}.1.part2.cadd.vcf.gz \
vcfs/by_chrom/${sample}.2.part1.cadd.vcf.gz \
vcfs/by_chrom/${sample}.2.part2.cadd.vcf.gz \
vcfs/by_chrom/${sample}.3.part1.cadd.vcf.gz \
vcfs/by_chrom/${sample}.3.part2.cadd.vcf.gz \
vcfs/by_chrom/${sample}.4.part1.cadd.vcf.gz \
vcfs/by_chrom/${sample}.4.part2.cadd.vcf.gz \
vcfs/by_chrom/${sample}.5.part1.cadd.vcf.gz \
vcfs/by_chrom/${sample}.5.part2.cadd.vcf.gz \
vcfs/by_chrom/${sample}.6.cadd.vcf.gz \
vcfs/by_chrom/${sample}.7.cadd.vcf.gz \
vcfs/by_chrom/${sample}.8.cadd.vcf.gz \
vcfs/by_chrom/${sample}.9.cadd.vcf.gz \
vcfs/by_chrom/${sample}.10.cadd.vcf.gz \
vcfs/by_chrom/${sample}.11.cadd.vcf.gz \
vcfs/by_chrom/${sample}.12.cadd.vcf.gz \
vcfs/by_chrom/${sample}.13.cadd.vcf.gz \
vcfs/by_chrom/${sample}.14.cadd.vcf.gz \
vcfs/by_chrom/${sample}.15.cadd.vcf.gz \
vcfs/by_chrom/${sample}.16.cadd.vcf.gz \
vcfs/by_chrom/${sample}.17.cadd.vcf.gz \
vcfs/by_chrom/${sample}.18.cadd.vcf.gz \
vcfs/by_chrom/${sample}.19.cadd.vcf.gz \
vcfs/by_chrom/${sample}.20.cadd.vcf.gz \
vcfs/by_chrom/${sample}.21.cadd.vcf.gz \
vcfs/by_chrom/${sample}.22.cadd.vcf.gz \
vcfs/by_chrom/${sample}.X.cadd.vcf.gz \
vcfs/by_chrom/${sample}.Y.cadd.vcf.gz | bgzip > vcfs/CADD/${sample}.vcf.gz


tabix -p vcf vcfs/CADD/${sample}.vcf.gz

rm vcfs/by_chrom/${sample}*


}
#=======================================
#END LARGE FUNCTION


# loop through samples
for i in `seq 1 $num_samples`
do

	sample=`echo $all_samples | cut -f $i -d' '`
	echo `date`
	echo "--working on ${i}/${num_samples} $sample"

	vcf=vcfs/gnomad/${sample}.vcf.gz

	if [ -e vcfs/CADD/$sample.vcf.gz.tbi ]
	then
		# do nothing
		echo 'already done'
	else
		annotate_cadd $sample
	fi

done

echo `date` fin

























