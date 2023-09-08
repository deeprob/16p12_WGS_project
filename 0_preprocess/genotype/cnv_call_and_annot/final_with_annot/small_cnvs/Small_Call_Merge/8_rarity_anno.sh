# Need to annotate DELs and DUPs separately

# DELs
# Annotate for intracohort frequency
bedtools intersect -a bed_files/7_dels.bed -b bed_files/7_dels.bed -f 0.5 -r -c > bed_files/8_dels.bed

# DUPs
# Intracohort frequency
bedtools intersect -a bed_files/7_dups.bed -b bed_files/7_dups.bed -f 0.5 -r -c > bed_files/8_dups.bed

# We will annotate gnomAD SV later as the annotation is a little more complicated
