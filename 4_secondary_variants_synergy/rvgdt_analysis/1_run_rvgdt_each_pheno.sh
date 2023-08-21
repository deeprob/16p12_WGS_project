#!/bin/bash

# originally prepared by jianyu
# edited by anastasia
# further edited by deepro
# edits
# 1. rvgdt cloned from github (git clone https://github.com/statgenetics/rv-gdt.git) and kept in directory, no need to edit PATH variable

## NOTE: This script is used to repeat RV-GDT runs several times to find stable significant genes across runs, so all files prepared for RV-GDT analysis can be found in previous analysis

# export PATH=/data5/16p12_WGS/TDT/rv-gdt/:$PATH # See edits 1

input_files_dir="/data5/deepro/wgs_16p/rvgdt/data/input_files"
intermediate_files_dir="/data5/deepro/wgs_16p/rvgdt/data/intermediate_files"
tmp_files_dir="/data5/deepro/wgs_16p/rvgdt/data/tmp_files"
results_dir="/data5/deepro/wgs_16p/rvgdt/data/results"

# input files which will be used in this script

# previously generated intermediate files which will be used in this script
genelist="${intermediate_files_dir}/rvgdt_genelist.txt"

# RV-GDT
for i in `seq 1 5`; do
mkdir -p "${results_dir}/${i}"
cd "${results_dir}/${i}"
for id in `cat "${input_files_dir}"/phenotypes2.id`; do
mkdir -p ./rvgdt_output_${id}
cd ./rvgdt_output_${id}
cat $genelist | parallel --jobs 64 rvgdt {}.${id} --geno $intermediate_files_dir/genos/{}.geno.txt --ped $intermediate_files_dir/pedigrees/family_pedigree.${id}.txt --max_iter 100
cd ..
done
cd ..
done


# merge results into a single file
ls -d */* | grep rvgdt_output | parallel "cat {}/* > {}.txt"
