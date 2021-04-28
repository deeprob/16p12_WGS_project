# CNV Postprocessing for RNAseq paper

For the RNAseq paper, we integrated the variants calls from the individual callers (CNVnator, Delly, Lumpy, Manta, and PennCNV) using custom scripts.

### parse_large_cnvs.ipynb

Formats large CNVnator calls (> 50 kbp) and PennCNV calls.

### parse_small_cnvs_by_caller.ipynb

Formats small CNV calls from CNVnator, Manta, Delly, and Lumpy.

### merge_small_cnv_caller.ipynb

Merges small CNV calls from CNVnator, Manta, Delly, and Lumpy.

### parse_small_cnvs.ipynb

Parses small CNVs and adds annotatios.

### combine_large_and_small_cnvs.ipynb

Combines large and small CNVs into one file.




