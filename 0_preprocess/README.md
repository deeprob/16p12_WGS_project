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

TODO: add from where the annotation datasets were obtained.

## CNV call and annotation
### WGS
Scripts used to call structural variants from WGS data divided by software.

#### CNVnator
*cnvnator.sh*: CNVnator calling script.

#### Delly
*0_delly.sh*:  Delly calling script.

*1_delly_merge.sh*: Merge calls into a single file. 

*2_delly_regenotype.sh*: Genotype from merged file.

*3_delly_merge.sh*: Merged samples into a single vcf file.

*4_delly_filter.sh*: Apply germline filter for delly.

*5_filter.sh*: Filter low quality variant calls.

#### Manta
*manta.sh*: Manta calling script.

#### Smoove lumpy
*1_smoove.sh*: Smoove calling script.

*2_union.sh*: Merge calls into a single vcf.

*3_regenotype.sh*: Genotype from merged file.

*4_get_failed_job_ids.sh*: Internal use for failed jobs.

*5_redo_failed_jobs.sh*: Internal use for failed jobs.

*6_paste.sh*: paste all the single sample VCFs with the same number of variants to get a single, squared, joint-called file.

*7_filter.sh*: Filter low quality variants. 

### Microarray
Scripts used to call structural variants from microarray data.

TODO: Find PennCNV script.

### Final with annotation
Finalize calls based on information from multiple callers and annotate them.

*1_combine_all_calls.py*: Large and Small CNV calls combined to a single file.

#### Large CNVs
Large CNV calls were a union of PennCNV and CNVnator calls.

##### PennCNV
