#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=batch4_gatk_pipe_array
#SBATCH -o ./logs/%a_gatk_pipe.log
#SBATCH -e ./logs/%a_gatk_error.log
#SBATCH --array=80
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=240:00:00
#SBATCH --mem-per-cpu=25G
#SBATCH --workdir /data4/16p12_WGS/batch4_1019

echo "Job started on `hostname` at `date`"

#Get list of all files
bulk_list=batch4_samples.txt
filename=`cat $bulk_list | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var'`

cd $filename

export TRIMMOMATIC=/data/software/Trimmomatic-0.36
export BWADIR=/data/software/bwa-0.7.13/
#export SAMDIR=/data/software/samtools-1.8
export GATKDIR=/data/software/gatk/gatk-3.7/
export PICARD=/data/software/picard-tools-2.9.0
export REF=/data/bx_references/hg19/ucsc.hg19.fasta
export DBSNP=/data/bx_references/hg19/dbsnp_138.hg19.vcf
export OMNI=/data/bx_references/hg19/1000G_omni2.5.hg19.vcf
export HAPMAP=/data/bx_references/hg19/hapmap_3.3.hg19.vcf
export MILLS=/data/bx_references/hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf
export G1K=/data/bx_references/hg19/1000G_phase1.indels.hg19.sites.vcf
export FWD=${filename}_R1.fastq.gz
export REV=${filename}_R2.fastq.gz

java -jar $TRIMMOMATIC/trimmomatic-0.36.jar PE -threads 10 -phred33 $FWD $REV ${filename}_R1_trim.fq.gz ${filename}_R1_unpaired.fq.gz ${filename}_R2_trim.fq.gz ${filename}_R2_unpaired.fq.gz LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20

$BWADIR/bwa mem -t 10 -R "@RG\tID:${filename}\tSM:${filename}\tPL:illumina\tLB:lib1" $REF ${filename}_R1_trim.fq.gz ${filename}_R2_trim.fq.gz | samtools sort -@ 10 -m 20G -o ${filename}_sorted.bam -
samtools view -@ 10 -h -m 24G -b -T $REF -o ${filename}_out.bam ${filename}_out.sam
samtools sort -@ 10 -m 20G -o ${filename}_sorted.bam ${filename}_out.bam
samtools index ${filename}_sorted.bam
rm ${filename}_out.sam
rm ${filename}_out.bam

java -Xmx240g -Djava.io.tmpdir=tmp -jar $PICARD/picard.jar MarkDuplicates \
I=${filename}_sorted.bam \
O=${filename}_rmdup.bam \
REMOVE_DUPLICATES=true \
M=${filename}_metrics.txt

samtools index ${filename}_rmdup.bam

java -Xmx240g -jar $GATKDIR/GenomeAnalysisTK.jar \
 -T BaseRecalibrator \
 -R $REF \
 -nct 10 \
 -I ${filename}_rmdup.bam \
 -knownSites $DBSNP \
 -knownSites $MILLS \
 -knownSites $G1K \
 -o ${filename}_BQSR_report.table

java -Xmx240g -jar $GATKDIR/GenomeAnalysisTK.jar \
 -T PrintReads \
 -R $REF \
 -nct 10 \
 -I ${filename}_rmdup.bam \
 -BQSR ${filename}_BQSR_report.table \
 -o ${filename}_rmdup_BQSR.bam

java -Xmx240g -jar $GATKDIR/GenomeAnalysisTK.jar \
 -T HaplotypeCaller \
 -R $REF \
 --dbsnp $DBSNP \
 -I ${filename}_rmdup_BQSR.bam \
 -o ${filename}_gvcf.g.vcf \
 --genotyping_mode DISCOVERY \
 -mbq 20 \
 --emitRefConfidence GVCF

echo "job completed at `date`"
