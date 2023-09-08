# Reduce columns for bed input
cut -f 1-3,6,9,19 call_tables/3_annotate_gencode.txt | tail -n+2 > tmp.bed

# Use bedtools to find all CNVs with a 50% reciprocal overlap of a given CNV for use in the python script
bedtools intersect -a tmp.bed -b tmp.bed -f 0.5 -r -loj > bed_files/4_inheritance_lookup.bed

rm tmp.bed
