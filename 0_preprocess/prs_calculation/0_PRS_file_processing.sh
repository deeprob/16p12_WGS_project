#Adapted from: https://choishingwan.github.io/PRS-Tutorial/base/

#Stage 1: Quality control of GWAS summary statistics files
#Filter for imputation INFO score >0.8 (when available); skip MAF>0.01 filter, as most GWAS statistics don't have minor allele frequency
gunzip -c sumstats_neuro_sum_ctg_format.txt.gz |awk 'NR==1 || ($14 > 0.8) {print}' |gzip  > filtering/neuroticism_info_filter.txt.gz

#Filter out ambiguous SNPs (SNPs where the alternate allele is complementary to the reference allele)
gunzip -c neuroticism_info_filter.txt.gz | awk '!( ($3=="A" && $4=="T") || ($3=="T" && $4=="A") || ($3=="G" && $4=="C") || ($3=="C" && $4=="G")) {print}' | gzip > neuroticism_ambig_snp_filter.txt.gz

#Remove duplicate SNPs (if any) from datasets
gunzip -c neuroticism_ambig_snp_filter.txt.gz |awk '{ print $2}' |sort | uniq -d > dup_snp.txt
gunzip -c neuroticism_ambig_snp_filter.txt.gz |grep -vf dup_snp.txt |gzip - > ../final_files/neuroticism_final.txt.gz

#Stage 2: Generation and QC of PLINK genotype files from microarray data

#Convert Illumina final report (with SNP ID, position, chromosome, and alleles) to lgen file
#Rscripts modified from: https://gist.github.com/RyanSchu/301ea0a77a21414391b54193dfcea9e0
Rscript illumina_to_lgen.R --illumina ../omniexpress_2019Dec/girirajan_grc_omni_3_191210_FinalReport.txt --outputdir ./

#Create "dummy" family file with internal sample IDs and sex (no family or phenotype info)
Rscript sex_to_fam.R

#Rename output files
mv mets.fam omniexpress_2019Dec.fam
mv mets.lgen omniexpress_2019Dec.lgen
mv mets.map omniexpress_2019Dec.map

#Convert LGEN/PED/MAP filesets to BED/BIM/FAM filesets for each batch
/data5/software/plink --lfile 16p12_hudson --out 16p12_hudson
/data5/software/plink --lfile 16p12_uw --out 16p12_uw
/data5/software/plink --lfile 16p12_yale --out 16p12_yale
/data5/software/plink --lfile omniexpress_2019 --out omniexpress_2019
/data5/software/plink --lfile omniexpress_2019Dec --out omniexpress_2019Dec

#Merge BED files together to create one fileset with all samples
/data5/software/plink --merge-list merge_filesets.txt --out 16p12_microarray_merged

#Flip strands on HudsonAlpha and Yale microarray data to match UW microarray data
/data5/software/plink --bfile 16p12_hudson --flip 16p12_microarray_merged.missnp --make-bed --out 16p12_hudson_flip
/data5/software/plink --bfile 16p12_yale --flip 16p12_microarray_merged.missnp --make-bed --out 16p12_yale_flip
/data5/software/plink --merge-list merge_filesets.txt --out 16p12_microarray_merged_flip

#Remove excluded samples, and change internal family and sample IDs to PSU/SG codes
/data5/software/plink --bfile 16p12_microarray_merged_flip --keep 16p12_microarray_samples.keep.txt --make-bed --out 16p12_microarray_filter
/data5/software/plink --bfile 16p12_microarray_filter --update-ids 16p12_microarray_samples.update.txt --make-bed --out 16p12_microarray_final

#Change sex information for select samples
/data5/software/plink2 --bfile 16p12_microarray_final --update-sex sample_sex_info.txt --make-bed --out 16p12_microarray_final_sex

#PLINK QC filters: Cohort MAF <0.05, Hardy-Weinberg equilibrium <1e-6 (under selection), geno<0.01 (SNPs missing in >1% of subjects), mind<0.01 (Individuals with >1% missing genotypes)
/data5/software/plink --bfile 16p12_microarray_final_sex --maf 0.05 --hwe 1e-6 --geno 0.01 --make-bed --out 16p12_microarray_final_QC
/data5/software/plink --bfile 16p12_microarray_final_QC --mind 0.01 --write-snplist --make-bed --out 16p12_microarray_final_FAM_QC
#8 samples removed at this stage: SG002, SG014, SG015, SG023, SG025, SG046, SG311, SG318
#Except for SG311, all samples were also filtered out of PennCNV for poor QC metrics

#Remove samples with extremely high heterozygosity
#Select SNPs in 200kb/50 variant windows that have LD r^2 scores of >0.25
/data5/software/plink --bfile 16p12_microarray_final_FAM_QC --keep 16p12_microarray_final_FAM_QC.fam --extract 16p12_microarray_final_FAM_QC.snplist --indep-pairwise 200 50 0.25 --out 16p12_microarray_final_FAM_QC
#Calculate F coefficient samples to estimate heterozygosity rate in each sample
/data5/software/plink --bfile 16p12_microarray_final_FAM_QC --extract 16p12_microarray_final_FAM_QC.prune.in --keep 16p12_microarray_final_FAM_QC.fam --het --out 16p12_microarray_final_FAM_QC
#Remove samples with >3SD heterozygosity rate, using the following R script:
	library(data.table)
	dat <- fread("16p12_microarray_final_FAM_QC.het")
	# Get samples with F coefficient within 3 SD of the population mean
	valid <- dat[F<=mean(F)+3*sd(F) & F>=mean(F)-3*sd(F)] 
	invalid <- dat[F>mean(F)+3*sd(F) | F<mean(F)-3*sd(F)]
	# print FID and IID for valid samples
	fwrite(valid[,c("FID","IID")], "16p12_microarray_final_FAM_QC.valid.sample", sep="\t") 
	fwrite(invalid[,c("FID","IID","F")], "16p12_microarray_final_FAM_QC.invalid.sample", sep="\t")
#7 samples were removed: SG163, SG165* (PSU056); SG170, SG172*, SG174* (PSU055); SG182, SG183* (PSU061) *=high
#No overlap with WGS samples

#Check sex of samples: Generate sexcheck file in PLINK
/data5/software/plink --bfile 16p12_microarray_final_FAM_QC --extract 16p12_microarray_final_FAM_QC.prune.in --keep 16p12_microarray_final_FAM_QC.valid.sample --check-sex --out 16p12_microarray_final_FAM_QC
#All samples have matching sexes--so skip the removal step
#Skip step to remove closely-related individuals

#Calculate principal component eigenvectors for population stratification
/data5/software/plink --bfile 16p12_microarray_final_FAM_QC --extract 16p12_microarray_final_FAM_QC.prune.in --pca 6 --out 16p12_microarray_final_FAM_QC

#Remove duplicate SNPs from PLINK files
cut -f 2 16p12_microarray_final_FAM_QC.bim | sort | uniq -d > 16p12_dup_snps.txt
/data5/software/plink --bfile 16p12_microarray_final_FAM_QC --exclude 16p12_dup_snps.txt --make-bed  --write-snplist --out 16p12_microarray_dedup

#For imputation step--create two sets of data, one with all QC-passed data and one with European samples only
#See separate .sh script for file processing


#Copy BIM file to retain original for downstream strand-flipping analysis
cp 16p12_microarray_dedup.bim 16p12_microarray_dedup_original.bim

#Adjust mismatched SNPs using strand-flipping in R
#Note: do this separately for each GWAS summary statistic dataset tested
library(data.table)
library(magrittr)
bim <- fread("16p12_microarray_dedup_original.bim") %>%
    setnames(., colnames(.), c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")) %>%
    #Change alleles to upper cases
    .[,c("B.A1","B.A2"):=list(toupper(B.A1), toupper(B.A2))]
# Read in GWAS summary statistic data (require data.table v1.12.0+)
gwas <- fread("../GWAS_statistics/final_files/Tourettes_final.txt.gz", fill=T) %>%
    #setnames(., colnames(.), c("CHR", "SNP", "A1", "A2", "BP", "INFO", "OR", "SE", "P", "NGT")) %>% 
    .[,c("CHR"):=list(strtoi(CHR))] #%>% #Include if column needs recast
    #.[,c("A1","A2"):=list(toupper(A1), toupper(A2))]
qc <- fread("16p12_microarray_dedup.snplist", header=F, fill=T)
#Merge data and filter for QC SNPs
info <- merge(bim, gwas, by=c("SNP", "CHR", "BP")) %>%
    .[SNP %in% qc[,V1]]
# Function for finding the complementary allele
complement <- function(x) {switch (x, "A" = "T", "C" = "G", "T" = "A", "G" = "C", return(NA)) }
# Get SNPs that have the same alleles across base and target
info.match <- info[A1 == B.A1 & A2 == B.A2, SNP]
# Identify SNPs that are complementary between base and target, and update BIM file
com.snps <- info[sapply(B.A1, complement) == A1 & sapply(B.A2, complement) == A2, SNP]
bim[SNP %in% com.snps, c("B.A1", "B.A2") := list(sapply(B.A1, complement), sapply(B.A2, complement))]
# Identify SNPs that need recoding and update BIM file
recode.snps <- info[B.A1==A2 & B.A2==A1, SNP]
bim[SNP %in% recode.snps, c("B.A1", "B.A2") := list(B.A2, B.A1)]
# identify SNPs that need recoding & complement, and update BIM file
com.recode <- info[sapply(B.A1, complement) == A2 & sapply(B.A2, complement) == A1, SNP]
bim[SNP %in% com.recode, c("B.A1", "B.A2") :=list(sapply(B.A2, complement), sapply(B.A1, complement))]
# Write the updated bim file
fwrite(bim, "16p12_microarray_dedup.bim", col.names=F, sep="\t")
#Output SNPs that do not match with target file
mismatch <- bim[!(SNP %in% info.match | SNP %in% com.snps | SNP %in% recode.snps | SNP %in% com.recode), SNP]
write.table(mismatch, "16p12_microarray_Tourettes_mismatch.txt", quote=F, row.names=F, col.names=F)

#Output final file
/data5/software/plink --bfile 16p12_microarray_dedup --make-bed --keep 16p12_microarray_final_FAM_QC.valid.sample --out 16p12_microarray_Tourettes_final --extract 16p12_microarray_dedup.snplist

#Stage 3: Calculate PRS using PLINK Clump method

#Log-transform ORs of summary statistics in R (skip this step for continuous traits)
	library(data.table)
	dat <- fread("schizophrenia_final.txt", fill=T)
	fwrite(dat[,OR:=log(OR)], "schizophrenia_final_transformed.txt", sep="\t")

#Perform clump step to cluster SNPs that are in linkage disequilibrium and are associated with the phenotype of interest
#Parameters: p1=select all SNPs for clumping, r2=remove SNPs with r2 index (haplotype MLE) >0.1, clump SNPs within 250kbp of target SNP
/data5/software/plink --bfile 16p12_WGS_ADHD_final --clump-p1 1 --clump-r2 0.1 --clump-kb 250 --clump adhd_final_transformed.txt.gz --clump-snp-field SNP --clump-field P --out 16p12_ADHD
#Extract index SNPs selected for downstream analysis
awk 'NR!=1{print $3}' 16p12_ADHD.clumped > 16p12_ADHD.valid.snp

#Generate file for ranges of p-values to test (do this once)
echo "0.001 0 0.001" > range_list.txt
echo "0.05 0 0.05" >> range_list.txt
echo "0.1 0 0.1" >> range_list.txt
echo "0.2 0 0.2" >> range_list.txt
echo "0.3 0 0.3" >> range_list.txt
echo "0.4 0 0.4" >> range_list.txt
echo "0.5 0 0.5" >> range_list.txt

#Generate input SNP list file for PRS calc (SNP ID and p-value)
zcat adhd_final_transformed.txt.gz | awk '{print $2,$9}' > adhd.SNP.pvalue.txt

#Perform PRS calculation
#Denote SNP ID, effective allele (A1), and effect size of GWAS statistic in --score parameter
/data5/software/plink --bfile 16p12_WGS_ADHD_final --score adhd_final_transformed.txt.gz 2 4 6 header --q-score-range range_list.txt adhd.SNP.pvalue.txt --extract 16p12_ADHD.valid.snp --out 16p12_WGS_ADHD
