#!/bin/bash


# concat all tables
> vcfs/formatted/ssc_d1.txt
cat samples_d1.txt | while read sample
do
	infile=vcfs/quality_filtered/${sample}.txt
	cat $infile >> vcfs/formatted/ssc_d1.txt
done


# format and select columns
python3 format_and_select_columns.py

# get strict inheritence
out=vcfs/formatted/ssc_d1.select_annotations.with_strict_inheritence.txt
python3 get_inheritence_strict.py > $out



