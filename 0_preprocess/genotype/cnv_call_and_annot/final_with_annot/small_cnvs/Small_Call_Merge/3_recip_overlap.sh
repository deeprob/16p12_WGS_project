#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=small_cnvs
#SBATCH -o logs/3_recip_overlap/batch_%a.log
#SBATCH -e logs/3_recip_overlap/batch_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=2G
#SBATCH --chdir /data5/16p12_WGS/structural_variants/sv_caller_postprocessing/2022_03_02/small_cnv_merge/
#SBATCH --array 1-287%10

echo `date` starting job on $HOSTNAME

SAMPLE=`head -n $SLURM_ARRAY_TASK_ID small_cnv_samples.list | tail -n1`
echo $SAMPLE

# Make a folder for each sample
OUTDIR=recip_overlap_bed_files/$SAMPLE
echo $OUTDIR
mkdir $OUTDIR

# Merge calls across individuals using 50% reciprocal overlap
# Keep calls that are present in at least 2 callers

CNVNATOR=sample_bed_files/cnvnator/$SAMPLE'_cnvnator_cnvs.bed'
MANTA=sample_bed_files/manta/$SAMPLE'_manta_cnvs.bed'
DELLY=sample_bed_files/delly/$SAMPLE'_delly_cnvs.bed'
LUMPY=sample_bed_files/lumpy/$SAMPLE'_lumpy_cnvs.bed'

# Start by annotating the number of times each call has a 50% reciprocal overlap in the other call sets
# CNVnator-Manta
bedtools intersect -a $CNVNATOR -b $MANTA -wa -f 0.5 -r -c > $OUTDIR/cnvnator_manta.bed
# Add Delly
bedtools intersect -a $OUTDIR/cnvnator_manta.bed -b $DELLY -wa -f 0.5 -r -c > $OUTDIR/cnvnator_manta_delly.bed
rm $OUTDIR/cnvnator_manta.bed
# Add Lumpy
bedtools intersect -a $OUTDIR/cnvnator_manta_delly.bed -b $LUMPY -wa -f 0.5 -r -c > $OUTDIR/cnvnator_manta_delly_lumpy.bed
rm $OUTDIR/cnvnator_manta_delly.bed

# Manta-CNVnator
bedtools intersect -a $MANTA -b $CNVNATOR -wa -f 0.5 -r -c > $OUTDIR/manta_cnvnator.bed
# Add Delly
bedtools intersect -a $OUTDIR/manta_cnvnator.bed -b $DELLY -wa -f 0.5 -r -c > $OUTDIR/manta_cnvnator_delly.bed
rm $OUTDIR/manta_cnvnator.bed
# Add Lumpy
bedtools intersect -a $OUTDIR/manta_cnvnator_delly.bed -b $LUMPY -wa -f 0.5 -r -c > $OUTDIR/manta_cnvnator_delly_lumpy.bed
rm  $OUTDIR/manta_cnvnator_delly.bed

# Delly-CNVnator
bedtools intersect -a $DELLY -b $CNVNATOR -wa -f 0.5 -r -c > $OUTDIR/delly_cnvnator.bed
# Add Manta
bedtools intersect -a $OUTDIR/delly_cnvnator.bed -b $MANTA -wa -f 0.5 -r -c > $OUTDIR/delly_cnvnator_manta.bed
rm $OUTDIR/delly_cnvnator.bed
# Add Lumpy
bedtools intersect -a $OUTDIR/delly_cnvnator_manta.bed -b $LUMPY -wa -f 0.5 -r -c > $OUTDIR/delly_cnvnator_manta_lumpy.bed
rm $OUTDIR/delly_cnvnator_manta.bed

# Lumpy-CNVnator
bedtools intersect -a $LUMPY -b $CNVNATOR -wa -f 0.5 -r -c > $OUTDIR/lumpy_cnvnator.bed
# Add Manta
bedtools intersect -a $OUTDIR/lumpy_cnvnator.bed -b $MANTA -wa -f 0.5 -r -c > $OUTDIR/lumpy_cnvnator_manta.bed
rm $OUTDIR/lumpy_cnvnator.bed
# Add Delly
bedtools intersect -a $OUTDIR/lumpy_cnvnator_manta.bed -b $DELLY -wa -f 0.5 -r -c > $OUTDIR/lumpy_cnvnator_manta_delly.bed
rm $OUTDIR/lumpy_cnvnator_manta.bed

echo `date` finished
