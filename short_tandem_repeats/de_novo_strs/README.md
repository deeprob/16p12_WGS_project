# De Novo STR with Monstr

Getting de novo STR with calls from GangSTR 2.5

### create_one_ped_per_trio.py

Creates one ped file per trio.

### create_one_vcf_per_trio.sh

Creates one VCF file per trio using the GangSTR 2.5 output.

### monstr.sh

Runs MonSTR on 16p12 cohort.

### monstr_chrX.sh

Runs MonSTR on chromosome X.

### monstr_long.sh

Runs monSTR based on flanking reads. monstr_long_chrX.sh does the same on chromsome X.

### format_and_filter.sh

Formats, filters, and annotates monSTR output with the following filters:

```
posterior probability > 0.5
locations with number of de novo mutations > 5 std. dev.
locations homozygous for the de novo allele
removes two samples with poor quality/contaminated sequencing
```

and the helper files:

```
filter_de_novo_strs.py
filter_de_novo_strs_long.py
filter_homozygous.py
prep_annovar.py
append_annovar.py
```



