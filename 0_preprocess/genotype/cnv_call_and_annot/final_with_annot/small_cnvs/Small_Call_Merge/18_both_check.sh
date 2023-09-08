#!/bin/bash

grep 'both' bed_files/16_inheritance.bed > bed_files/18_both_cnvs.bed

lines=`wc -l bed_files/18_both_cnvs.bed | cut -f 1 -d ' '`

# Make graphs for 29 de novo plots
for (( line=1; line<=$lines; line++ ));
do
	# Get information from CNV call
	SAMPLE=`head -n $line bed_files/18_both_cnvs.bed | tail -n1 | cut -f 7`
	CHR=`head -n $line bed_files/18_both_cnvs.bed | tail -n1 | cut -f 1`
	START=`head -n $line bed_files/18_both_cnvs.bed | tail -n1 | cut -f 2`
	END=`head -n $line bed_files/18_both_cnvs.bed | tail -n1 | cut -f 3`
	TYPE=`head -n $line bed_files/18_both_cnvs.bed | tail -n1 | cut -f 4 | sed 's/<//' | sed 's/>//'`

	# Get parent information
	FATHER=`grep -P "PSU[0-9][0-9][0-9]\t$SAMPLE" ../../2022_02_14/Feb_2022.fam | cut -f 3`
	MOTHER=`grep -P "PSU[0-9][0-9][0-9]\t$SAMPLE" ../../2022_02_14/Feb_2022.fam | cut -f 4`

	# Get final calls for all family members as a BED file for plotting
	grep -e $SAMPLE -e $FATHER -e $MOTHER bed_files/16_inheritance.bed | cut -f 1-3,13 | tail -n+2 | bedtools sort | bgzip -c > bed_files/call_track.bed.gz
	tabix bed_files/call_track.bed.gz

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
		-o 'both_plots/'$SAMPLE'_'$CHR'_'$START'_'$END'_'$TYPE'.png' \
		-c $CHR \
		-s $START \
		-e $END \
		-t $TYPE \
		-A bed_files/call_track.bed.gz

	# Remove annotation track
	rm bed_files/call_track.bed.gz
done

echo `date` finished
