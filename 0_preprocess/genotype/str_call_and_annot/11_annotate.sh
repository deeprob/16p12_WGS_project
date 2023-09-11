#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=annovar
#SBATCH -o logs/1_annotate.log
#SBATCH -e logs/1_annotate.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=30G
#SBATCH --chdir /data5/16p12_WGS/structural_variants/expansions_strs2/
#SBATCH --nodelist ramona

# Prep to get statistics
# I ran this on Ramona


# run annovar
infile=/data5/16p12_WGS/structural_variants/expansions_strs/output/16p12_cohort.expansions_2SD.annovar_input.bed
outprefix=/data5/16p12_WGS/structural_variants/expansions_strs2/16p12_cohort.expansions_2SD

perl /data4/software/annovar/table_annovar.pl $infile  /data5/bx_reference/hg19/annotations/annovar \
 -buildver hg19 \
 -out $outprefix \
 -remove \
 -protocol refGene,wgEncodeGencodeBasicV19 \
 -operation g,g \
 -nastring . \
 -arg '-hgvs','-hgvs'





