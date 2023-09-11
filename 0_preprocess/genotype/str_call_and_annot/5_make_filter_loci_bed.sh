#!/bin/bash
# creates a bed file of loci with population level filters


#======================================================
# Autosomal chromosomes
#======================================================
vcf=/data5/16p12_WGS/structural_variants/vcf_callers/output/mergestr_by_family/mergestr_filtered.vcf.gz
out=output/merged_loci_remove.bed

vcf-query $vcf -f '%CHROM\t%POS\t%INFO/END\t%FILTER\n' |  awk '{if ($4 != "PASS") print $0 ;}' >  $out

#======================================================
# For Chromosome X there are two VCFs with filters
#======================================================

vcf=/data5/16p12_WGS/structural_variants/vcf_callers/output/mergestr_by_family_chrX/mergestr_filtered.vcf.gz
out=output/merged_loci_remove_chrX.bed

vcf-query $vcf -f '%CHROM\t%POS\t%INFO/END\t%FILTER\n' |  awk '{if ($4 != "PASS") print $0 ;}' >  $out

vcf=/data5/16p12_WGS/structural_variants/vcf_callers/output/mergestr_by_family_chrX/mergestr_female_only_filtered.vcf.gz
out=output/merged_loci_remove_chrX_hw.bed

vcf-query $vcf -f '%CHROM\t%POS\t%INFO/END\t%FILTER\n' |  awk '{if ($4 != "PASS") print $0 ;}' >  $out

# combine the two bed files
cat output/merged_loci_remove_chrX.bed output/merged_loci_remove_chrX_hw.bed | sort > output/merged_loci_remove_chrX_combined.bed

# check that the lines sum up
wc -l output/merged_loci_remove_chrX.bed 
wc -l output/merged_loci_remove_chrX_hw.bed 
wc -l output/merged_loci_remove_chrX_combined.bed 

