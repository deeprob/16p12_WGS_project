#!/bin/bash

ls peds | while read ped
do
	echo $ped
	child=`echo $ped | cut -f1 -d.`
	family=`head -n1 peds/$ped | cut -f1`
	samples=`cat peds/$ped | cut -f2`
	
	samples=`echo $samples | sed 's/\s\+/,/g'`
	
	echo $family
	echo $samples

	vcf=../gangstr/output/gangstr/${family}.vcf.gz

	bcftools view -s $samples $vcf | bgzip > vcfs/$child.vcf.gz
	tabix -p vcf vcfs/$child.vcf.gz

done

# Chromosome X ==============================
ls peds | while read ped
do
	echo $ped
	child=`echo $ped | cut -f1 -d.`
	family=`head -n1 peds/$ped | cut -f1`
	samples=`cat peds/$ped | cut -f2`
	
	samples=`echo $samples | sed 's/\s\+/,/g'`
	
	echo $family
	echo $samples

	vcf=../gangstr/output/gangstr/${family}.chrX.vcf.gz

	bcftools view -s $samples $vcf | bgzip > vcfs/$child.chrX.vcf.gz
	tabix -p vcf vcfs/$child.chrX.vcf.gz

done
# ==============================


echo 'done'
