# Need to annotate DELs and DUPs separately

# DELs
# First annotate samples with gnomAD SV count
#bedtools intersect -a bed_files/8_dels.bed -b /data5/16p12_WGS/structural_variants/sv_caller_postprocessing/2022_02_14/gnomad_v2.1_sv.dels.bed -wa -f 0.5 -r -c > bed_files/9.1_dels.bed
# Filtering for gnomAD SV frequency will be complicated, so we will do that in another step

# Then, annotate for microarray count
bedtools intersect -a bed_files/8_dels.bed -b /data5/16p12_WGS/structural_variants/sv_caller_postprocessing/2022_02_14/cnv_controls_hg19_dels.txt -f 0.5 -r -c > bed_files/9.2_dels.bed

# Finally, annotate for intracohort frequency
bedtools intersect -a bed_files/9.2_dels.bed -b bed_files/9.2_dels.bed -f 0.5 -r -c > bed_files/9.3_dels.bed

# DUPs
# gnomAD SV
#bedtools intersect -a bed_files/8_dups.bed -b /data5/16p12_WGS/structural_variants/sv_caller_postprocessing/2022_02_14/gnomad_v2.1_sv.dups.bed -wa -f 0.5 -r -c > bed_files/9.1_dups.bed
# Microarray controls
bedtools intersect -a bed_files/8_dups.bed -b /data5/16p12_WGS/structural_variants/sv_caller_postprocessing/2022_02_14/cnv_controls_hg19_dups.txt -f 0.5 -r -c > bed_files/9.2_dups.bed
# Intracohort frequency
bedtools intersect -a bed_files/9.2_dups.bed -b bed_files/9.2_dups.bed -f 0.5 -r -c > bed_files/9.3_dups.bed
