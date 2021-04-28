#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=spark_d1_inh
#SBATCH -o logs/strict_inheritence.log
#SBATCH -e logs/strict_inheritence.err
#SBATCH --cpus-per-task=2
#SBATCH --time=240:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --chdir /data5/SPARK_WES/annotations


echo `date` started on $HOSTNAME

out=vcfs/formatted/ssc_d1.select_annotations.with_strict_inheritence.txt
python3 get_inheritence_strict.py > $out






echo `date` finished





