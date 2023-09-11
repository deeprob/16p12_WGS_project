library(data.table)
library(magrittr)
bim <- fread("16p12_microarray_final_hg19_original.bim") %>%
    setnames(., colnames(.), c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2"))
qc <- fread("16p12_microarray_final_hg19.snplist", header=F, fill=T)

# Read in GWAS summary statistic data (require data.table v1.12.0+)
#Template:
gwas <- fread("./final_files/Tourettes_final.txt.gz", fill=T) %>%
    #setnames(., colnames(.), c("CHR", "SNP", "A1", "A2", "BP", "INFO", "OR", "SE", "P", "NGT")) %>% 
    .[,c("CHR"):=list(strtoi(CHR))] #%>% #Include if column needs recast
    #.[,c("A1","A2"):=list(toupper(A1), toupper(A2))]

#Separate commands for each GWAS summary statistic file:

gwas <- fread("./final_files/educational_attainment_final.txt.gz", fill=T) %>%
    setnames(., colnames(.), c("SNP","CHR", "BP", "A1", "A2", "EAF", "BETA", "SE", "P")) %>% 
    .[,c("CHR"):=list(strtoi(CHR))]

gwas <- fread("./final_files/intelligence_final.txt.gz", fill=T) %>%
    setnames(., colnames(.), c("SNP", "UNIQ_ID", "CHR", "BP", "A1", "A2", "EAF", "Z", "BETA", "SE", "P", "N", "INFO", "EFFECT_DIR")) %>% 
    .[,c("CHR"):=list(strtoi(CHR))] %>%
    .[,c("A1","A2"):=list(toupper(A1), toupper(A2))]

gwas <- fread("./final_files/schizophrenia_final.txt.gz", fill=T) %>%
    setnames(., colnames(.), c("CHR", "SNP", "A1", "A2", "BP", "INFO", "OR", "SE", "P", "NGT")) %>% 
     .[,c("CHR"):=list(sub('...', '', CHR))] %>% .[,c("CHR"):=list(strtoi(CHR))] #Remove "chr" from text file

gwas <- fread("./final_files/autism_final.txt.gz", fill=T) %>%
    .[,c("CHR"):=list(strtoi(CHR))]

gwas <- fread("./final_files/cross_disorder_final.txt.gz", fill=T) %>%
    setnames(., colnames(.), c("CHR", "BP", "SNP", "A2", "A1", "BETA", "SE", "P", "NGT", "FCAS", "FCON", "IMPINFO","NEFFDIV2", "NCAS", "NCON", "DIRE","ASSET","m.SCZ","m.BIP","m.MD","m.ASD","m.ADHD","m.TS","m.AN","m.OCD")) %>% 
    .[,c("CHR"):=list(strtoi(CHR))]

#Merge data and filter for QC SNPs
#NOTE: did not merge by chr/bp, as base pairs of SNPs did not agree w/ each other (diff. HG builds?)
info <- merge(bim, gwas, by=c("SNP")) %>% 
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
fwrite(bim, "16p12_microarray_final_hg19.bim", col.names=F, sep="\t")
#Output SNPs that do not match with target file
mismatch <- bim[!(SNP %in% info.match | SNP %in% com.snps | SNP %in% recode.snps | SNP %in% com.recode), SNP]
write.table(mismatch, "16p12_cross_disorder_mismatch.txt", quote=F, row.names=F, col.names=F)

quit()

#Unix command following R script
/data5/software/plink --bfile 16p12_microarray_final_hg19 --make-bed --out 16p12_microarray_final_cross_disorder --extract 16p12_microarray_final_hg19.snplist



#Commands for other PRS datasets
gwas <- fread("./GWAS_final_files/adhd_final.txt.gz", fill=T) %>%
    .[,c("CHR"):=list(strtoi(CHR))]

gwas <- fread("./GWAS_final_files/alcohol_dependence_final.txt.gz", fill=T) %>%
    .[,c("CHR"):=list(strtoi(CHR))]

gwas <- fread("./GWAS_final_files/bipolar_final.txt.gz", fill=T) %>%
    .[,c("CHR"):=list(strtoi(CHR))]

gwas <- fread("./GWAS_final_files/BMI_final.txt.gz", fill=T) %>%
    setnames(., colnames(.), c("CHR", "BP", "A2", "A1", "SNP", "GMAF", "MAF", "BETA", "SE", "P")) %>% 
    .[,c("CHR"):=list(strtoi(CHR))]

gwas <- fread("./GWAS_final_files/height_final.txt.gz", fill=T) %>%
    setnames(., colnames(.), c("CHR", "BP", "A2", "A1", "SNP", "GMAF", "MAF", "BETA", "SE", "P")) %>% 
    .[,c("CHR"):=list(strtoi(CHR))]

gwas <- fread("./GWAS_final_files/insomnia_final.txt.gz", fill=T) %>%
    .[,c("CHR"):=list(strtoi(CHR))]


gwas <- fread("./GWAS_final_files/MDD_final.txt.gz", fill=T) %>%
    .[,c("CHR"):=list(strtoi(CHR))]

gwas <- fread("./GWAS_final_files/neuroticism_final.txt.gz", fill=T) %>%
    setnames(., colnames(.), c("SNP", "ALT_ID", "CHR", "BP", "A1", "A2", "EAF", "MAF", "Z", "BETA", "SE", "P", "N", "INFO")) %>% 
    .[,c("CHR"):=list(strtoi(CHR))]

gwas <- fread("./GWAS_final_files/OCD_final.txt.gz", fill=T) %>%
    .[,c("CHR"):=list(strtoi(CHR))]

gwas <- fread("./GWAS_final_files/PTSD_final.txt.gz", fill=T) %>%
    .[,c("CHR"):=list(strtoi(CHR))]


gwas <- fread("./GWAS_final_files/Tourettes_final.txt.gz", fill=T) %>%
    .[,c("CHR"):=list(strtoi(CHR))]