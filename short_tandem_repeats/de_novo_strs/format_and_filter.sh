#!/bin/bash

# concat all files
ls output/monstr/ | grep locus_summary | grep -v chrX | grep -v all.all | while read l; do cat output/monstr/$l; done > output/monstr/all.locus_summary.tab
ls output/monstr/ | grep locus_summary | grep chrX | grep -v all.all | while read l; do cat output/monstr/$l; done > output/monstr/all.chrX.locus_summary.tab
ls output/monstr/ | grep all_mutations | grep -v chrX | grep -v all.all | while read l; do cat output/monstr/$l; done > output/monstr/all.all_mutations.tab
ls output/monstr/ | grep all_mutations | grep chrX | grep -v all.all | while read l; do cat output/monstr/$l; done > output/monstr/all.chrX.all_mutations.tab

# long STR calls
ls output/monstr/long/ | grep locus_summary | grep -v chrX | grep -v all.all | while read l; do cat output/monstr/long/$l; done > output/monstr/long/all.locus_summary.tab
ls output/monstr/long/ | grep locus_summary | grep chrX | grep -v all.all | while read l; do cat output/monstr/long/$l; done > output/monstr/long/all.chrX.locus_summary.tab
ls output/monstr/long/ | grep all_mutations | grep -v chrX | grep -v all.all | while read l; do cat output/monstr/long/$l; done > output/monstr/long/all.all_mutations.tab
ls output/monstr/long/ | grep all_mutations | grep chrX | grep -v all.all | while read l; do cat output/monstr/long/$l; done > output/monstr/long/all.chrX.all_mutations.tab



#=================================================
# filter for de novo strs
# these scripts remove headers that are in the middle of the text file
# and only keeps de novo mutations with posterior probability > 0.5
python3 filter_de_novo_strs.py output/monstr/all.all_mutations.tab > output/filtered/all.all_mutations.filtered.tab
python3 filter_de_novo_strs.py output/monstr/all.chrX.all_mutations.tab > output/filtered/all.chrX.all_mutations.filtered.tab

# long STRs
python3 filter_de_novo_strs_long.py output/monstr/long/all.all_mutations.tab > output/filtered/all.all_mutations.long.filtered.tab
python3 filter_de_novo_strs_long.py output/monstr/long/all.chrX.all_mutations.tab > output/filtered/all.chrX.all_mutations.long.filtered.tab


#=================================================

# remove sex chromosomes from all mutations
cat output/filtered/all.all_mutations.filtered.tab | awk '$1 != "chrX" && $1 != "chrY"' > output/filtered/all.all_mutations.filtered.no_sex.tab

# append autosomal and X chromosome
cat output/filtered/all.all_mutations.filtered.no_sex.tab > output/filtered/all.filtered.tab
tail -n +2 output/filtered/all.chrX.all_mutations.filtered.tab >> output/filtered/all.filtered.tab

# long STRs 
# long STRsdon't have any sex chromosomes, so just move them to new file
cat output/filtered/all.all_mutations.long.filtered.tab > output/filtered/all.long.filtered.tab


#=================================================

# filter out SG473 who has > 1000 mutations
# and filter out SG136 whose mother SG138 WGS sample is bad
cat output/filtered/all.filtered.tab | grep -v SG473 | grep -v SG136 > output/filtered/all.filtered2.tab


# more filters
# filters out children with number of de novo mutations > 5 std dev.
# filters out locations with number of de novo mutations > 5 std dev.
python3 /data5/software/STRDenovoTools-1.0.0/scripts/qc_denovos.py \
--all-mutations-file output/filtered/all.filtered2.tab \
--filtered-mutations-file output/filtered/all.filtered3.tab \
--log-file logs/qc_de_novo.log \
--filter-loc-denovos 5

# filter if the child is homozygous for the de novo allele
python3 filter_homozygous.py output/filtered/all.filtered3.tab ../gangstr/all_batches_fixed.ped >output/filtered/all.filtered4.tab

# other filters if both children have the same de novo and contraction and expansion bias I won't do

# annovar annotations
# all STR annotations are here:
# /data5/16p12_RNA/scripts/str/strs.full_anno.hg19_multianno.txt 
python3 prep_annovar.py output/filtered/all.filtered4.tab output/filtered/annovar_input.tsv


# run ANNOVAR
perl /data4/software/annovar/table_annovar.pl output/filtered/annovar_input.tsv  /data4/software/annovar/humandb/ \
 -buildver hg19 \
 -out output/filtered/annovar \
 -remove \
 -protocol refGene \
 -operation gx \
 -nastring . \
 -thread 6 \
 -xref /data4/software/annovar/humandb/gene_annotations.txt \
 -arg '-hgvs'



# add annovar annotations to original file
outfilename=output/filtered/all.annotated.tab
python3 append_annovar.py output/filtered/all.filtered4.tab output/filtered/annovar.hg19_multianno.txt $outfilename




