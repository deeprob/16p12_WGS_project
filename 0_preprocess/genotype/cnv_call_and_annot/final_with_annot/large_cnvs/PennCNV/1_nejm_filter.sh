# Convert PennCNV calls to BED file
cut -f 1-3,6,9 -d , ../../2022_02_14/PennCNV_calls_combined_final.csv | sed 's/,/\t/g' | tail -n+2 > tmp.bed

# use bedtools to find regions that have a 50% reciprocal overlap with NEJM CNVs
bedtools intersect -wa -wb -a tmp.bed -b ../../2022_02_14/nejm_cnv_hg19.bed -f 0.5 -r > bed_files/1_nejm_filter.bed

rm tmp.bed
