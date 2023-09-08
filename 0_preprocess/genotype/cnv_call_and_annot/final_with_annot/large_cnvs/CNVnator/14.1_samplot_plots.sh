#!/bin/bash

# Make graphs for every large CNV

line=$1

# Get information from CNV call
SAMPLE=`head -n $line bed_files/13_frequency_brandler_filter.bed | tail -n1 | cut -f 7`
CHR=`head -n $line bed_files/13_frequency_brandler_filter.bed | tail -n1 | cut -f 1`
START=`head -n $line bed_files/13_frequency_brandler_filter.bed | tail -n1 | cut -f 2`
END=`head -n $line bed_files/13_frequency_brandler_filter.bed | tail -n1 | cut -f 3`
TYPE=`head -n $line bed_files/13_frequency_brandler_filter.bed | tail -n1 | cut -f 4 | sed 's/<//' | sed 's/>//'`

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
