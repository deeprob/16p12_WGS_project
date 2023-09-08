# Annotate the gnomAD SV frequency

# Create a header line
gnomad_head=`head -1 /data5/16p12_WGS/structural_variants/sv_caller_postprocessing/2022_02_14/gnomad_v2.1_sv.sites.bed`
echo  -e 'Chr\tStart\tEnd\tType\tName\tLength\tSample\tIntracohort_count\t'$gnomad_head > 10_header.txt
sed -i 's/ \+/\t/g' 10_header.txt

cat 10_header.txt > bed_files/10_gnomadSV_anno_dels.bed
cat 10_header.txt > bed_files/10_gnomadSV_anno_dups.bed

rm 10_header.txt

# Because of the way the gnomAD SVs are labelled, we'll have to do a more complicated annotation than the one used for intracohort
# This step builds a lookup table of CNVs and all gnomAD CNVs it is associated with
bedtools intersect -a bed_files/9_dels_intracohort_filter.bed -b /data5/16p12_WGS/structural_variants/sv_caller_postprocessing/2022_02_14/gnomad_v2.1_sv.dels.bed -loj -f 0.5 -r >> bed_files/10_gnomadSV_anno_dels.bed

bedtools intersect -a bed_files/9_dups_intracohort_filter.bed -b /data5/16p12_WGS/structural_variants/sv_caller_postprocessing/2022_02_14/gnomad_v2.1_sv.dups.bed -loj -f 0.5 -r >> bed_files/10_gnomadSV_anno_dups.bed

