# Describes all the files and their locations required to run rvgdt

## Input files (required)
1. Participants info file

    Filename: ./input_files/16p12_All_Participants_v9.xlsx

    Original Location: Dropbox/16p12.2 project/Human patients project/WGS paper/3_Cohort information/16p12_All_Participants_v9.xlsx

2. Structural Variant file

    Filename: ./input_files/sv_calls_combined.txt

    Original Location: Dropbox/16p12.2 project/Human patients project/WGS paper/7_Structural variants/sv_calls_combined.txt

3. Child phenotypic domain information

    Filename: ./input_files/prelim_child_domains_Jul16.xlsx 

    Original Location: Dropbox/16p12.2 project/Human patients project/Phenotypic analysis/de Vries & Psychiatric/WGS/MC_Children/prelim_child_domains_Jul16.xlsx

4. ?

    Filename: ./input_files/Phenotype_splits_new.xlsx 

    Original Location: /data5/anastasia/rvgdt/Phenotype_splits_new.xlsx

5. ?

    Filename: ./input_files/phenotypes2.id

    Original Location: /data5/anastasia/rvgdt/phenotypes2.id

6. Pedigree file (ext change)

    Filename: ./input_files/all_batches.ped

    Original Location: Dropbox/16p12.2 project/Human patients project/Pedigree Information/Mar_2022.fam 

7. Single nucleotide rare variant file 

    Filename: ./input_files/16p12_TDT_missense_rare_filter.vcf

    Original Location: /data5/16p12_WGS/annotations/pvcf_annotations/bvcfs/16p12_TDT_missense_rare_filter.vcf

8. Genome fasta file

    Filename: /data/bx_references/hg19/ucsc.hg19.fasta

    Original Location: /data/bx_references/hg19/ucsc.hg19.fasta


## Intermediate files (generated)
1. Parsed participants file

    Filename: ./intermediate_files/tdt_sample_list.txt

2. Sample list with trios file

    Filename: ./intermediate_files/tdt_sample_list_trios.txt

3. Known pedigree file

    Filename: ./intermediate_files/all_batches.no_unkowns.ped

4. Merged and filtered vcf file

    Filename: ./intermediate_files/16p12_WGS_genotype_coding_snv_rare_filterSV.vcf

5. left-normalized file

    Filename: ./intermediate_files/16p12_WGS_genotype_coding_snv_rare_filterSV_norm.vcf

6. GT only vcf

    Filename: ./intermediate_files/16p12_WGS_genotype_coding_snv_rare_filterSV.GT.vcf

7. rvgdt gene list

    Filename: ./intermediate_files/rvgdt_genelist.txt
