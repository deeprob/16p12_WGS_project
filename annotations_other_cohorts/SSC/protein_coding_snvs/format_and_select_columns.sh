cat format_and_select_columns.sh 
#!/bin/bash


# format and select annotations


# concat annotations for all chromosomes
for i in `seq 1 22` X Y
do
	cat vcfs/protein_coding/chr${i}.txt
done > vcfs/protein_coding/all_chromosomes.txt


# select annotations
python3 format_and_select_columns.py 

# add LEOUF
python3 add_leouf.py


