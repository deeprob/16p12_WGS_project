#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=vcf2table
#SBATCH -o logs/8_vcf2table/batch_%a.log
#SBATCH -e logs/8_vcf2table/batch_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=10G
#SBATCH --chdir /data5/16p12_WGS/annotations/2022_02_04/
#SBATCH --array 2-346
#SBATCH --nodelist qingyu

echo `date` starting job on $HOSTNAME


export PATH=/data5/anastasia/sw/bcftools-1.12:$PATH


SAMPLE=`cat samples.csv | head -n $SLURM_ARRAY_TASK_ID | tail -n1 | cut -f2 -d,`
INPUT_VCF=/data5/16p12_WGS/annotations/2022_02_04/exonic_vcfs/cadd/$SAMPLE.vcf.gz
OUTPUT_FILE=/data5/16p12_WGS/annotations/2022_02_04/exonic_variants/by_sample/$SAMPLE.txt


echo $SAMPLE
echo input file: $INPUT_VCF
echo output file: $OUTPUT_FILE




bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%Func.refGene\t%Gene.refGene\t%GeneDetail.refGene\t%ExonicFunc.refGene\t%AAChange.refGene\t%Func.wgEncodeGencodeBasicV19\t%Gene.wgEncodeGencodeBasicV19\t%GeneDetail.wgEncodeGencodeBasicV19\t%ExonicFunc.wgEncodeGencodeBasicV19\t%AAChange.wgEncodeGencodeBasicV19\t%gnomad_exome_AF\t%gnomad_genome_AF\t%CADD_PHRED\t%CADD_RawScore[\t%SAMPLE\t%GT\t%DP\t%AD\t%GQ\t%PL\t%RGQ\t%SB]\n' $INPUT_VCF > $OUTPUT_FILE
# bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%Func.refGene\t%Gene.refGene\t%GeneDetail.refGene\t%ExonicFunc.refGene\t%AAChange.refGene\t%gnomad_exome_AF\t%gnomad_genome_AF\t%CADD_PHRED\t%CADD_RawScore[\t%SAMPLE\t%GT\t%DP\t%AD\t%GQ\t%PL\t%RGQ\t%SB]\n' $INPUT_VCF > $OUTPUT_FILE






echo `date` finished


