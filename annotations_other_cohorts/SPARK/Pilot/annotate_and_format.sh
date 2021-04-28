#!/bin/bash


export PATH=/data5/anastasia/sw/bcftools-1.12/build/bin/:$PATH

# left norm
in=/data4/SPARK/Pilot/VCF/SPARK_pilot1379ind.ColumbiaJointCall.vcf.gz
fasta=/data/bx_references/GRCh37/human_g1k_v37.fasta
out=SPARK_pilot1379ind.ColumbiaJointCall.norm.vcf.gz
bcftools norm -f $fasta -m-both $in | bgzip --threads 4 > $out
tabix -p vcf $out



in=SPARK_pilot1379ind.ColumbiaJointCall.norm.vcf.gz
outprefix=SPARK_pilot1379ind.ColumbiaJointCall.annovar
perl /data4/software/annovar/table_annovar.pl $in /data4/software/annovar/humandb/ \
 -buildver hg19 \
 -out $outprefix \
 -remove \
 -protocol refGene,generic \
 -operation gx,f \
 -nastring . \
 -thread 6 \
 -vcfinput \
 -xref  /data5/SPARK_WES/annotations/release_june2020/gene_annotations.txt \
 -genericdbfile mpc_values_v2_final.txt \
 -arg '-hgvs',

bgzip $outprefix.hg19_multianno.vcf
tabix -p vcf $outprefix.hg19_multianno.vcf.gz


# add MPC
/data5/software/slivar  expr --vcf SPARK_pilot1379ind.ColumbiaJointCall.annovar.hg19_multianno.vcf.gz \
	-g /data5/anastasia/liftover/hg19/mpc_hg19_gnotate.zip \
	| bgzip \
	> SPARK_pilot1379ind.ColumbiaJointCall.annovar.hg19_multianno.with_MPC.vcf.gz

# filter by locus and sample
in=SPARK_pilot1379ind.ColumbiaJointCall.annovar.hg19_multianno.with_MPC.vcf.gz
out=SPARK_pilot1379ind.ColumbiaJointCall.filtered.txt
python3 filter_by_locus_and_sample.py $in > $out

# format and select columns
python3 format_and_select_columns.py 

# strict inheritence
out=SPARK_pilot1379ind.ColumbiaJointCall.filtered.select_annotations.with_strict_inheritence.txt
python3 get_inheritence_strict.py > $out




