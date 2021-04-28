#!/bin/bash

# script to get, sort and index GRCh37 SegDups file

# Note: GRCh37GenomicSuperDup.sorted.gz was obtained from UCSC Table Browser
# with the options: Mammal, HUman, GRCh37/hg19, Repeats, Segemntal Dups, genomicSuperDups, output format BED
# https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1047622799_NYRO2yV9XWpcQdLYnWOEDqBIAVUJ&clade=mammal&org=&db=hg19&hgta_group=rep&hgta_track=genomicSuperDups&hgta_table=genomicSuperDups&hgta_regionType=genome&position=&hgta_outputType=primaryTable&hgta_outFileName=


wget "https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1047622799_NYRO2yV9XWpcQdLYnWOEDqBIAVUJ&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_genomicSuperDups&hgta_ctDesc=table+browser+query+on+genomicSuperDups&hgta_ctVis=pack&hgta_ctUrl=&fbQual=whole&fbUpBases=200&fbDownBases=200&hgta_doGetBed=get+BED"
mv 'hgTables?hgsid=1047622799_NYRO2yV9XWpcQdLYnWOEDqBIAVUJ&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_genomicSuperDups&hgta_ctDesc=table+browser+query+on+genomicSuperDups&hgta_ctVis=pack&hgta_ctUrl=&fbQual=whole&fbUpBases=200&fb' GRCh37genomicSuperDups.bed

bgzip GRCh37genomicSuperDups.bed 

bedtools sort -i GRCh37genomicSuperDups.bed.gz | bgzip > GRCh37genomicSuperDups.sorted.bed.gz

tabix -p vcf GRCh37genomicSuperDups.sorted.bed.gz 

rm GRCh37genomicSuperDups.bed.gz 
