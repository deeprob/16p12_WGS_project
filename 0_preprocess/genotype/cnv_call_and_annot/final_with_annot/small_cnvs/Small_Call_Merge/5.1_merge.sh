#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=small_cnvs_merge
#SBATCH -o logs/5_merge/batch_%a.log
#SBATCH -e logs/5_merge/batch_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=2G
#SBATCH --chdir /data5/16p12_WGS/structural_variants/sv_caller_postprocessing/2022_03_02/small_cnv_merge/
#SBATCH --array 1-287%50

echo `date` starting job on $HOSTNAME

SAMPLE=`head -n $SLURM_ARRAY_TASK_ID small_cnv_samples.list | tail -n1`
echo $SAMPLE

python 5_merge_calls.py $SAMPLE

echo `date` finished
