# use bedtools to find regions that have a 50% reciprocal overlap with NEJM CNVs
bedtools intersect -wa -wb -a bed_files/4_merge_files.bed -b ../../nejm_cnv_hg19.bed -f 0.5 -r > bed_files/5_nejm_filter.bed
