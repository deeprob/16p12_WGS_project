# Iterate through individual files and save to one large table

# Start with file header
echo -e 'Chr\tStart\tEnd\tType\tID\tLength\tSample' > bed_files/5_merge_files.bed

for i in {2..346}
do
	SAMPLE=`cat lumpy_files.csv | head -n $i | tail -n1 | cut -f2 -d,`

	# Check if we are using sample for WGS, if not skip
	if grep -q $SAMPLE /data5/16p12_WGS/structural_variants/sv_caller_postprocessing/2022_02_14/16p12_WGS_Consent_Participants.csv
	then
		BED_FILE='bed_files/4_size_filter/'$SAMPLE'_size_filter.bed'
		tail -n+2 $BED_FILE >> bed_files/5_merge_files.bed
	else
		echo Skipping $SAMPLE - not in participants table
	fi
done
