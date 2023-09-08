#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=separate_by_sample
#SBATCH -o logs/1_separate_by_sample/delly.log
#SBATCH -e logs/1_separate_by_sample/delly.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=10G
#SBATCH --chdir /data5/16p12_WGS/structural_variants/sv_caller_postprocessing/2022_03_02/small_cnv_merge/

# Need to separate CNV calls by sample for merging
CNVNATOR=/data5/16p12_WGS/structural_variants/sv_caller_postprocessing/2022_02_14/cnvnator/small_cnv_processing/bed_files/7_str_brandler_filter.bed
MANTA=/data5/16p12_WGS/structural_variants/sv_caller_postprocessing/2022_02_14/manta/bed_files/8_str_brandler_filter.bed
DELLY=/data5/16p12_WGS/structural_variants/sv_caller_postprocessing/2022_02_14/delly/bed_files/8_str_brandler_filter.bed

caller=$DELLY

line_counts=`wc -l $caller | cut -f 1 -d ' '`
echo $line_counts
for ((line=2; line <= $line_counts; line++));
do
	sample=`head -n $line $caller | tail -1 | cut -f 7`

	echo $line $caller $sample

	# Save line to file for sample
	head -n $line $caller | tail -1 >> sample_bed_files/delly/$sample'_delly_cnvs.bed'

done
