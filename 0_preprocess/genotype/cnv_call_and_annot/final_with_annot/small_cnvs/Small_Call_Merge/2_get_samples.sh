# Get sample names from log files
ls -lthr sample_bed_files/* | sed 's/ \+/\t/g' | cut -f 9 | cut -f 1 -d '_' | sort | uniq | tail -n+3 > small_cnv_samples.list

# Restrict to only 16p12.1 deletion families that are not excluded from WGS
