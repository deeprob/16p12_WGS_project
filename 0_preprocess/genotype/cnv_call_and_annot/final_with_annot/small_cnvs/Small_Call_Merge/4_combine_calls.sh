# Combine calls from all 3 callers into one file for each person
for line in {1..287};
do
	SAMPLE=`head -n $line small_cnv_samples.list | tail -n1`
	echo $SAMPLE

	# Remove calls that are only present in one caller (counts for other 2 callers are 0)
	grep -vP '\t0\t0\t0' recip_overlap_bed_files/$SAMPLE/cnvnator_manta_delly_lumpy.bed > 4_combined_calls/$SAMPLE'_combined.bed'
	grep -vP '\t0\t0\t0' recip_overlap_bed_files/$SAMPLE/manta_cnvnator_delly_lumpy.bed >> 4_combined_calls/$SAMPLE'_combined.bed'
	grep -vP '\t0\t0\t0' recip_overlap_bed_files/$SAMPLE/delly_cnvnator_manta_lumpy.bed >> 4_combined_calls/$SAMPLE'_combined.bed'
	grep -vP '\t0\t0\t0' recip_overlap_bed_files/$SAMPLE/lumpy_cnvnator_manta_delly.bed >> 4_combined_calls/$SAMPLE'_combined.bed'

done
