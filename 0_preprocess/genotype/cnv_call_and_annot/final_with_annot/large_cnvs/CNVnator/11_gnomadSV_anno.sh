# Annotate the gnomAD SV frequency

# Create a header line
gnomad_head=`head -1 /data5/16p12_WGS/structural_variants/sv_caller_postprocessing/2022_02_14/gnomad_v2.1_sv.sites.bed`
echo  -e 'Chr\tStart\tEnd\tType\tName\tLength\tSample\tMicroarray_count\tIntracohort_count\tmicroarray_freq\t'$gnomad_head > 11_header.txt
sed -i 's/ \+/\t/g' 11_header.txt

cat 11_header.txt > bed_files/11_gnomadSV_anno_dels.bed
cat 11_header.txt > bed_files/11_gnomadSV_anno_dups.bed

rm 11_header.txt

# Because of the way the gnomAD SVs are labelled, we'll have to do a more complicated annotation than the one used for microarray controls and intracohort
# This step builds a lookup table of CNVs and all gnomAD CNVs it is associated with
bedtools intersect -a bed_files/10_dels_frequency_filter.bed -b /data5/16p12_WGS/structural_variants/sv_caller_postprocessing/2022_02_14/gnomad_v2.1_sv.dels.bed -loj -f 0.5 -r >> bed_files/11_gnomadSV_anno_dels.bed

bedtools intersect -a bed_files/10_dups_frequency_filter.bed -b /data5/16p12_WGS/structural_variants/sv_caller_postprocessing/2022_02_14/gnomad_v2.1_sv.dups.bed -loj -f 0.5 -r >> bed_files/11_gnomadSV_anno_dups.bed

