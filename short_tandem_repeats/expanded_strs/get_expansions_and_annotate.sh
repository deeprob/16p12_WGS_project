#!/bin/bash

# Prep to get statistics
# I ran this on Ramona

# merge

vcfs=`cat ../gangstr/families_list.txt | sed -s 's/^/\/data5\/16p12_WGS\/structural_variants\/gangstr\/output\/vcfs_filtered\//g' | sed -s 's/$/.autosomal_and_sex.vcf.gz/g'`
# replace space with comma
vcfs=`echo $vcfs | sed 's\ \,\g'`


mergeSTR --vcfs $vcfs --out output/16p12_cohort

# vcf to table 
python3 vcf_to_table.py > output/16p12_cohort.tsv

# get stats

statSTR --vcf output/16p12_cohort.vcf.gz \
	--thresh \
	--afreq \
	--acount \
	--mean \
	--mode \
	--var \
	--numcalled \
	--out output/16p12_cohort.statstr

# get STR expansions > 2SD from mean
python3 get_expansions_2sd.py > output/16p12_cohort.expansions_2SD.tsv

# prep for annovar
tail +2 output/16p12_cohort.expansions_2SD.tsv | cut -f1-5 | uniq > output/16p12_cohort.expansions_2SD.annovar_input.bed

# run annovar
infile=output/16p12_cohort.expansions_2SD.annovar_input.bed
outprefix=output/16p12_cohort.expansions_2SD

perl /data4/software/annovar/table_annovar.pl $infile  /data4/software/annovar/humandb/ \
 -buildver hg19 \
 -out $outprefix \
 -remove \
 -protocol refGene,knownGene \
 -operation gx,gx, \
 -nastring . \
 -thread 6 \
 -xref /data4/software/annovar/humandb/gene_annotations.txt \
 -arg '-hgvs','-hgvs'


# combine annovar annotations with other table
python3 add_annotations.py

