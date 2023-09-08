# Because lumpy files are bigger than other callers, they were not condensed into one file
# Just copy each file into a new directory in this folder
for i in {2..346};
do
	SAMPLE=`head -n $i ../lumpy/lumpy_files.csv | tail -n 1 | cut -f 2 -d ,`
	echo $SAMPLE `date`

	# Check if we are using sample for WGS, if not skip
	if grep -q $SAMPLE /data5/16p12_WGS/structural_variants/sv_caller_postprocessing/2022_02_14/16p12_WGS_Consent_Participants.csv
	then
		cp ../lumpy/bed_files/8_str_overlap/$SAMPLE'_str_brandler_filter.bed' sample_bed_files/lumpy/$SAMPLE'_lumpy_cnvs.bed'
	else
		echo Skipping $SAMPLE - not in participants table
	fi
done
