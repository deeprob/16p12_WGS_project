# Use bedtools to find all CNVs with a 50% reciprocal overlap of a given CNV for use in the python script
bedtools intersect -a bed_files/16_all_calls_combined.bed -b bed_files/16_all_calls_combined.bed -f 0.5 -r -loj > bed_files/17.1_inheritance_lookup.bed

