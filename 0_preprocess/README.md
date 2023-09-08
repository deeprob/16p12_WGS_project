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
*1_nejm_filter.sh*: Find regions that have a 50% reciprocal overlap with NEJM CNVs

*2_filter_penncnv.py*: Filter PennCNV calls based on frequency, gene number and SegDups and Centromere/Telomere overlap

*3_annotate_gencode.py*: Annotate gencode genes

*4_inheritance_lookup.sh*: Find all CNVs with a 50% reciprocal overlap of a given CNV

*5_add_inheritance.py*: Add inheritance information for CNVs

*6_filter_genic.py*: Filter calls for only those that overlap genes

*7_confirm_denovo.sh*: Visually confirm the de novo calls using the visualizations Matt had made

*8_update_inheritance.py*: Update inheritance after manual verification

*9_finalize_calls.py*: Limit calls to only those we are using for WGS analysis making sure all samples passed microarray

*10_separate_genes.py*: Separate calls into gene-level calls

*11_annotate_genes.py*: Get LOEUF and OMIM annotations for genes

##### CNVnator
*1_get_filenames.sh*: Get all CNVnator called files. 

*2_adjacent_filter.sh*: Merge adjacent calls per individual for each caller

*3_split_by_size.sh*: Separate small and large CNV calls

*4_merge_files.sh*: Merge files into one big file. 

*5_nejm_filter.sh*: Find regions that have a 50% reciprocal overlap with NEJM CNVs

*6_low_confidence_region_filter.sh*: Find regions that have less than a 50% reciprocal overlap with centromeres, SegDups, regions of low mappability, and V(D)J recombination regions from Brandler et al. Science 2016

*7_str_filter.sh*: Remove any CNV with a breakpoint in an STR region

*8_type_split.sh*: Split CNVs into Dels and Dups for easy frequency filtering

*9_rarity_anno.sh*: Annotate dels and dups with population and cohort frequency

*10_frequency_filter.py*: remove samples with microarray controls frequency > 0.1 or intracohort frequency > 10

*10.1_frequency_anno.py*: Annotate with microarray frequency

*11_gnomadSV_anno.sh*: Annotate the gnomAD SV frequency

*12_gnomadSV_filter.py*: Filter the calls using the gnomAD SV AF

*12.1_gnomadSV_anno.py*: Filter the calls using the gnomAD SV AF

*13_combine_dels_dups.sh*: Combine deletion and duplication calls

*13.1_combine_nejm_dels_dups.sh*: Combine deletion and duplication calls

*14_samplot_plots.sh*: Make graphs for every large CNV

*14.1_samplot_plots.sh*: Make graphs for every large CNV

*14.2_nejm_samplot_plots.sh*: Make graphs for every large CNV

*15_annotate_gencode.py*: Annotate gencode genes using gencode.v19.parsed.exons.csv

*16_combine_calls.py*: ombined NEJM and other calls into one callset

*17_add_inheritance.py*: Add cnv inheritance info

*17.1_inheritance_lookup.sh*: Find all CNVs with a 50% reciprocal overlap of a given CNV

*18_denovo_check.sh*: Visually confirm the de novo calls using the visualizations Matt had made

*19_microarray_cnvnator_overlap.sh*: Get the overlap of the microarray calls and the CNVnator calls

*19.1_microarray_cnvnator_overlap.py*: Get the overlap of the microarray calls and the CNVnator calls

*20_finalize_calls.py*: Finalize CNV calls based on overlap.

*21_separate_genes.py*: Separate calls into gene-level calls

*22_annotate_genes.py*: Get LOEUF and OMIM annotations for genes


#### Small CNVs
Small CNV calls were called by at least two of CNVnator, Delly, Lumpy, Manta.

##### CNVnator


##### Delly

##### Lumpy

##### Manta

##### Small Call Merge