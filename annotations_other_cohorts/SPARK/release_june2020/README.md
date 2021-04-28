# Annotating SSC Rare Variants

### filter_ref.sh

Removes variant calls that are ref or missing to reduce file size

### left_norm.sh

Left normalizes variant calls and trims ALT alleles.

### annovar.sh

Anntates with annovar.

### annotate_gnomad.sh

Annotates gnomad population frequency.

### annotate_cadd.sh

Annotates CADD score.

### annotate_MPC.sh

Annotates MPC.

### filter_quality.sh

Filters VCF for quality and a few other metrics. 
```
QUAL >= 50
gnomad freq <= 0.001
location in exonic region or splice site
variants function is splic, lof, missense, or frameshift
per-sample quality metrics
```

### format_and_select_columns.sh

Concats all variant calls into one file, formats into a table, and selects columns

### get_inheritence_strict.sh

Gets strict inheritence. For strict inheritence, a sample's inheritence is left as missing unless both parents pass per-sample quality metrics for the loci.

