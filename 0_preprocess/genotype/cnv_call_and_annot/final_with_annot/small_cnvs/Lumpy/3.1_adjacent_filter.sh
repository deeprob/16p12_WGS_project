#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=adjacent_filter
#SBATCH -o logs/3_adjacent_filter/batch_%a.log
#SBATCH -e logs/3_adjacent_filter/batch_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=2G
#SBATCH --chdir /data5/16p12_WGS/structural_variants/sv_caller_postprocessing/2022_03_02/lumpy
#SBATCH --array 2-346%50

echo `date` starting job on $HOSTNAME

SAMPLE=`cat lumpy_files.csv | head -n $SLURM_ARRAY_TASK_ID | tail -n1 | cut -f2 -d,`
INPUT_FILE='vcfs/2_get_cnvs/'$SAMPLE'_dels_dups.vcf'

echo $SAMPLE
echo input file: $INPUT_FILE

python 3_adjacent_filter.py $INPUT_FILE $SAMPLE

echo `date` finished

# NOTES:
# Job 273 failed (disk space/memory error?)
