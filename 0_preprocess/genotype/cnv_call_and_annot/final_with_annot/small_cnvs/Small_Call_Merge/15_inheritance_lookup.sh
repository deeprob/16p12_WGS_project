# Make an inheritance lookup table - similar to the gnomAD SV table
bedtools intersect -a bed_files/14_gene_filter.bed -b bed_files/14_gene_filter.bed -f 0.5 -r -loj > bed_files/15_inheritance_lookup.bed
