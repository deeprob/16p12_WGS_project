# STR variant calling
GangSTR (variant calling), dumpSTR (variant filtering), and mergeSTR (VCF merging and statistics) for 16p12 cohort.

### gangstr_by_family.sh

GangSTR variant calling for 16p12 deletion cohort (except family PSU008).

### gangstr_by_family_PSU008.sh

GangSTR variant calling for family PSU008.

### dumpstr_by_family.sh

Filters gangSTR variant calls with --min-call-DP 20	--max-call-DP 1000	--filter-spanbound-only	--filter-badCI	--readlen 150

Note: these are the same parameters as Gymrek's paper (excluding --require-support 2 which results in empty VCF)

--require-support 2 uses the AD field, which gangSTR doesn't output.

see https://www.nature.com/articles/s41586-020-03078-7#Sec6

### mergestr_by_family.sh

Merges VCFs output by dumpstr_by_family.sh with default MergeSTR parameters.

### get_GRCh37SuperDups.sh

Get the GRCh37genomicSuperDups.sorted.bed.gz file from the UCSC Table Browser.

### dumpstr_by_family2.sh

Run dumpSTR on merged VCF with the options --min-locus-hwep 0.00001 --min-locus-callrate 0.8 --filter-regions GRCh37genomicSuperDups.sorted.bed.gz	--filter-regions-names SEGDUP

## Calling STRs on Chromosome X

The following scripts are related to calling and filtering short tandem repeats on chromosome X.

### make_sex_file.py

Python script to create a file listing each samples ploidy for chromosome X. It takes in a pedigree file and outputs a file with the format

```
sample  ploidy
SG123 2
SG012 1
SG234 2
...
```

### gangstr_by_family_chromX.sh

GangSTR variant calling for 16p12 deletion cohort for chromosome X (except family PSU008). This script is the same as gangstr_by_family.sh except it uses the option --ploidy.

### gangstr_by_family_chromX_PSU008.sh

GangSTR variant calling for chromosome X for family PSU008.

### dumpstr_by_family_chromX.sh

Filters gangSTR variant calls with --min-call-DP 10	--max-call-DP 1000	--filter-spanbound-only	--filter-badCI	--readlen 150

### mergestr_by_family_chromX.sh

Merges VCFs output by dumpstr_by_family_chromX.sh with default MergeSTR parameters.

### dumpstr_by_family2_chrX.sh

Run dumpSTR on merged VCF with the options --min-locus-callrate 0.8 --filter-regions GRCh37genomicSuperDups.sorted.bed.gz	--filter-regions-names SEGDUP

### dumpstr_by_family2_chrX_hw.sh

Run dumpSTR on merged VCF with the options --min-locus-hwep 0.00001 only on female samples.
