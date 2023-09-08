# use bedtools to find regions that have less than a 50% reciprocal overlap with centromeres, SegDups, regions of low mappability, and V(D)J recombination regions from Brandler et al. Science 2016
bedtools intersect -wa -a bed_files/5_merge_files.bed -b ../filtered_regions_hg19.bed -f 0.5 -r -v -header > bed_files/7_brandler_filter.bed
