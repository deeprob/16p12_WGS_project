# Get the overlap of the microarray calls and the CNVnator calls

# Make a BED file from the PennCNV calls
cut -f 1-3,6,9 -d , ../../PennCNV_calls_combined_final.csv | sed 's/,/\t/g' | tail -n+2 > bed_files/19_penncnv.bed

# Use bedtools for overlap
bedtools intersect -a bed_files/17_inheritance.bed -b bed_files/19_penncnv.bed -wa -wb -loj -f 0.5 -r >  bed_files/19_overlap.bed
