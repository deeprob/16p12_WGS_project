# Combine deletion and duplication calls
cat bed_files/12.1_nejm_dels.bed > bed_files/13_nejm.bed
tail -n+2 bed_files/12.1_nejm_dups.bed >> bed_files/13_nejm.bed

