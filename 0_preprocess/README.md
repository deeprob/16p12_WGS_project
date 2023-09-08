# Description
This directory holds the preprocessing scripts for raw data analyzed in the manuscript. It is further sub-divided into folders.

# Phenotype
*0_Child_Domains_Jul16.ipynb*: Preparation of different phenotypic domains based on observed phenotypes

# Genotype
## Genome alignment and SNV calls
*0_gatk_pipe_array.sh*: Trimming, aligning, filtering, duplicate removal and SNV calling using gatk's best practices pipeline.

## SNV annotation
*1_list_of_files.py*: Get all gvcf files created through gatk pipeline that needs to be annotated. 

*2_filter_quality.sh*: Filter variants for base quality score and read depth.

*3_annovar_genes.sh*: Get annovar gene annotations.

*4_filter_exonic.sh*: Filter to keep only exonic variants. 

*5_annotate_gnomad.sh*: Get gnomad annotations using vcf anno.

*6_filter_pop_freq.sh*: Filter by gnomad frequency.

*7_annotate_cadd.sh*: Get cadd annotations.

*8_vcf2table.sh*: Convert the vcf file to a table.

*9_combine_tables.sh*: Combine individual tables into a single table.

*10_filter_variant_types.py*: Filter by CADD Phred scores and mutation type. 

*11_intracohort_filter.py*: Filter variants that are present in less than 10 unique individuals in our cohort.

*12_gencode_gene_ids.py*: Map the genes to their symbols and ids. 

*13_louef_scores.py*: Assign loeuf scores to the genes.

*14_omim_annotations.py*: Get OMIM annotations for genes. 

*15_finalize_variants.sh*: Copy the final file to another file.
