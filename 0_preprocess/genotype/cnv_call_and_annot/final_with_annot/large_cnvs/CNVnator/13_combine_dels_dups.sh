# Combine deletion and duplication calls
cat bed_files/12_gnomadSV_filter_dels.bed > bed_files/13_frequency_brandler_filter.bed
tail -n+2 bed_files/12_gnomadSV_filter_dups.bed >> bed_files/13_frequency_brandler_filter.bed

