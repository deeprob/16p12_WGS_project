#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=samplot_plots
#SBATCH -o logs/12_samplot_plots/batch_%a.log
#SBATCH -e logs/12_samplot_plots/batch_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=2G
#SBATCH --chdir /data5/16p12_WGS/structural_variants/sv_caller_postprocessing/2022_02_14/cnvnator/large_cnv_processing/
#SBATCH --array 2-3

# Make graphs for every large CNV

echo `date` starting job on $HOSTNAME

export PATH=/afs/bx.psu.edu/user/c/ces6136/miniconda3/bin/samplot:$PATH

# Get information from CNV call
SAMPLE=`head -n $SLURM_ARRAY_TASK_ID bed_files/13_frequency_brandler_filter.bed | tail -n1 | cut -f 7`
CHR=`head -n $SLURM_ARRAY_TASK_ID bed_files/13_frequency_brandler_filter.bed | tail -n1 | cut -f 1`
START=`head -n $SLURM_ARRAY_TASK_ID bed_files/13_frequency_brandler_filter.bed | tail -n1 | cut -f 2`
END=`head -n $SLURM_ARRAY_TASK_ID bed_files/13_frequency_brandler_filter.bed | tail -n1 | cut -f 3`
TYPE=`head -n $SLURM_ARRAY_TASK_ID bed_files/13_frequency_brandler_filter.bed | tail -n1 | cut -f 4 | sed 's/<//' | sed 's/>//'`

# Get parent information
FATHER=`grep -P "PSU[0-9][0-9][0-9]\t$SAMPLE" ../../Feb_2022.fam | cut -f 3`
MOTHER=`grep -P "PSU[0-9][0-9][0-9]\t$SAMPLE" ../../Feb_2022.fam | cut -f 4`

# Get BAM files
SAMPLE_BAM=`grep $SAMPLE bam.list`
if [ "$FATHER" = 0 ]
then
	FATHER_BAM=''
else
	FATHER_BAM=`grep $FATHER bam.list`
fi
if [ "$MOTHER" = 0 ]
then
	MOTHER_BAM=''
else
	MOTHER_BAM=`grep $MOTHER bam.list`
fi

echo $SAMPLE $CHR $START $END $TYPE $FATHER $MOTHER
echo $SAMPLE_BAM $FATHER_BAM $MOTHER_BAM

if [ "$FATHER_BAM" = '' ]
then
	FATHER=''
fi
if [ "$MOTHER_BAM" = '' ]
then
	MOTHER=''
fi

echo $SAMPLE $CHR $START $END $TYPE $FATHER $MOTHER

# Make samplot plots
samplot plot \
	-n $SAMPLE $FATHER $MOTHER \
	-b $SAMPLE_BAM $FATHER_BAM $MOTHER_BAM \
	-o 'plots/'$SAMPLE'_'$CHR'_'$START'_'$END'_'$TYPE'.png' \
	-c $CHR \
	-s $START \
	-e $END \
	-t $TYPE

echo `date` finished
