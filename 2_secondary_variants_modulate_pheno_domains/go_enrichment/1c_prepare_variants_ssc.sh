#!/bin/bash



SSC_DIR="C:/Users/corny/Dropbox/16p12.2 project/Human patients project/WGS paper/31_SSC analysis/Genotype-phenotype integration/1_variant_preparation/"






for cohort in dbd_tier1_snvs large_rare_deletions large_rare_duplications nejm_deletions nejm_duplications
do
	echo $cohort
	
	cat "${SSC_DIR}/samples/${cohort}.csv" | tail -n +2 | cut -f1 -d, | sort | uniq > cohorts/SSC_${cohort}.csv
	
	
	mkdir -p variants/SSC/$cohort
	cp "${SSC_DIR}/variants/$cohort/rare_deleterious_snvs.csv" variants/SSC/$cohort
	cp "${SSC_DIR}/variants/$cohort/strs.csv" variants/SSC/$cohort/strs_std_2.csv
	cp "${SSC_DIR}/variants/$cohort/deletions.csv" variants/SSC/$cohort/dels.csv
	cp "${SSC_DIR}/variants/$cohort/duplications.csv" variants/SSC/$cohort/dups.csv
	
done











