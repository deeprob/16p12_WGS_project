# Redoing STR variant calling with GangSTR 2.5 

We originally ran GangSTR 2.4 on the 16p12 cohort. However, to to get de novo STR calls, we need the ENCLREADS parameter that comes with GangSTR 2.5.

Here, we rerun GangSTR on the locations that weren't filtered out by dumpSTR in the original run.

### make_filter_loci_bed.sh

Creates a bed file of variants that should be filtered out due to population-level filters made by DumpSTR.

### filter_vcf.sh

Filter VCF file with bed file made by make_filter_loci_bed.sh

### get_important_calls.sh

For each family, makes a bedfile of locations to keep

### filter_vcf_chrX.sh

Combines the steps of filter_vcf.sh and get_important_calls.sh for chromosome X.

### gangstr*.sh

Runs GangSTR 2.5 on the 16p12 cohort with the important calls.
Note: needed to switch to all_batches_fixed.ped file because of incorrect sex labels.



