# Combine all the files into one for the cohort

echo -e 'Chr\tStart\tEnd\tType\tID\tLength\tSample' > bed_files/6_combine_samples.bed

for i in {1..287};
do
	SAMPLE=`head -n $i small_cnv_samples.list | tail -1`
	tail -n+2 merged_calls/$SAMPLE'_merged.bed' >> bed_files/6_combine_samples.bed
done
