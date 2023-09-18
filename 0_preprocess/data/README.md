# Annotation files description
This file lists sources from where the annotation files were downloaded.

1. [Clinvar](#clinvar)
2. [GENCODE](#gencode)
3. [Gene annotations](#gene-annotations)
4. [Genehancer annotations](#genehancer-annotations)
5. [Gnomad annotations](#gnomad-annotations)
6. [OMIM annotations](#omim-annotations)

# Clinvar
All files were downloaded from here: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/ on February 14th, 2022

# Gencode
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz

# Gene annotations
Werling_NatGen_2018_Supp
	- Supplementary Tablles 1-13 from Werling et al. NatGen 2018
	- Downloaded on March 9, 2022
	- Download link: https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0107-y/MediaObjects/41588_2018_107_MOESM3_ESM.zip
	- Original paper: https://www.nature.com/articles/s41588-018-0107-y
	- Folder contains Excel sheets for every Supplementary Table from the paper
	- We are interested in SuppTable06_GeneLists.xlsx

Full-LoF-Table-Data.csv
	- Geisinger Developmental Brain Disorder (DBD) Gene Database LoF table download
	- Downloaded on March 9, 2022
	- Website: https://dbd.geisingeradmi.org/
	- Download link: https://dbd.geisingeradmi.org/downloads/Full-LoF-Table-Data.csv
	- Original paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5333489/
	- File contains number of LOF cases in 6 neurodevelopmental disorders (ID, ASD, epilepsy, ADHD, schizophrenia, and bipolar disorder) and "Tier"s of pathogenicity - see below
	- Tier meanings:
		- Tier 1: Genes with three or more de novo pathogenic loss-of-function variants (High confidence)
		- Tier 2: Genes with two de novo pathogenic loss-of-function variants (High confidence)
		- Tier 3: Genes with one de novo pathogenic loss-of-function variant (Emerging candidate)
		- Tier 4: Genes with no de novo pathogenic loss-of-function variant (Emerging candidate)
		- Tier AR: Genes with autosomal recessive inheritance (High confidence)

DDG2P.csv.gz
	- Gene2Phenotype DD gene-disease pairs and attributes
	- Downloaded on March 14, 2022
	- Website: https://www.ebi.ac.uk/gene2phenotype
	- Download link: https://www.ebi.ac.uk/gene2phenotype/downloads/DDG2P.csv.gz
	- Original paper: http://europepmc.org/abstract/MED/31147538
	- File contains high-confidence gene-disease pairs for developmental delay diseases, along with attrbites of potential disease-causing mutations
	
SFARI-Gene_genes_01-11-2022release_03-14-2022export.csv
	- SFARI Gene Database download
	- Downloaded on March 14, 2022
	- Website: https://gene.sfari.org/database/human-gene/
	- Download link: https://gene.sfari.org//wp-content/themes/sfari-gene/utilities/download-csv.php?api-endpoint=genes
	- Original paper: https://molecularautism.biomedcentral.com/articles/10.1186/2040-2392-4-36
	- File contains genes, autism gene scores, and other information
	- Gene score meanings: (see https://gene.sfari.org/about-gene-scoring/ for more information)
		- S: Syndromic
		- 1: High confidence
		- 2: Strong candidate
		- 3: Suggestive evidence
		
SZDB
	- SZDB2.0 GWAS genes, CNV genes, and Exome sequencing genes
	- Downloaded on March 16, 2022
	- Website: http://szdb.org/
	- Original paper: https://pubmed.ncbi.nlm.nih.gov/27451428/
	- Update paper: https://pubmed.ncbi.nlm.nih.gov/32385526/
	- File contains genes associated with schizophrenia from a variety of sources, see following links for details:
		- GWAS: http://szdb.org/gwas1.php
		- CNV: http://szdb.org/cnv_gene1.php
		- Exome: Data is a combination of these studies - http://szdb.org/exome-publication.php
		
Epilepsy gene list.xlsx
	- Genes associated with epilepsy from Wang Seizure 2017
	- This file was compiled by Matt(?) and represents all genes from all tables from that paper
	- Paper link: https://pubmed.ncbi.nlm.nih.gov/28007376/

RVIS_Unpublished_ExACv2_March2017.txt
	- RVIS scores
	- Downloaded on April 19, 2022
	- Website: https://genic-intolerance.org/
	- Download link: https://genic-intolerance.org/data/RVIS_Unpublished_ExACv2_March2017.txt
	- Original paper: http://www.plosgenetics.org/article/info:doi/10.1371/journal.pgen.1003709
	- File contains several columns of data pertaining to RVIS scores
		- The relevant column is %RVIS[pop_maf_0.05%(any)] (RVIS percentile)
	
# Genehancer annotations
genehancer_ucsc_download.txt
	- This file is all files in the Chromosome annotations folder combined

GeneHancer_Version_4-4.xlsx
	- This file was downloaded from GeneCards (https://www.genecards.org/cgi-bin/carddisp.pl?gene=IKZF1&keywords=genehancer)
		- Scroll down to GeneHancer section and click "Download" GeneHancer 2017 data
	- This data is from the original GeneHancer paper (https://academic.oup.com/database/article/doi/10.1093/database/bax028/3737828?login=false)

GeneaHancerDataSheet-2021
	- This file was downloaded from GeneCards
	- It describes the GeneHancer data
	
# Gnomad annotations
https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz

# Omim annotations
Files were downloaded from https://www.omim.org/downloads on January 31st, 2022
