#Code adapted from: https://choishingwan.github.io/PRS-Tutorial/ldpred/

library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

library(data.table)
library(magrittr)


#Step 1: Load phenotype + summary statistics files (not needed for auto model)
#Load combination phenotype + covariate file
#pheno <- fread('16p12_pheno_cov.txt')

#Load HapMap3 SNP dataset
info <- readRDS(url("https://github.com/privefl/bigsnpr/raw/master/data-raw/hm3_variants.rds"))
sumstats <- fread("../final_files/educational_attainment_final.txt.gz", fill=T) %>%
    setnames(., colnames(.), c("SNP","CHR", "BP", "A1", "A2", "EAF", "BETA", "SE", "P")) %>% 
    .[,c("CHR"):=list(strtoi(CHR))]

names(sumstats) <-c("rsid","chr","pos","a1","a0","EAF","beta","beta_se","p") #Change based on header
#Note: for quantitive traits, used total sample size as the effective sample size--not sure if this is correct
sumstats$n_eff <- 766345
sumstats <- sumstats[sumstats$rsid%in% info$rsid,] #Filter out SNPs not in HapMap3 dataset

#Step 2: Calculate LD matrix
# Get maximum amount of cores
NCORES <- nb_cores()
# Open a temporary file
tmp <- tempfile(tmpdir = "tmp-data")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
# Initialize variables for storing the LD score and LD matrix
corr <- NULL
ld <- NULL
# We want to know the ordering of samples in the bed file 
fam.order <- NULL

# preprocess the bed file (only need to do once for each data set)
#snp_readBed("../16p12_microarray_final_educational_attainment.bed")
# now attach the genotype object
obj.bigSNP <- snp_attach("../16p12_microarray_final_educational_attainment.rds")
# extract the SNP information from the genotype
map <- obj.bigSNP$map[-3]
names(map) <- c("chr", "snpid", "pos", "a1", "a0")
map <- subset(map,chr>0) #Remove SNPs w/o location information
map<-subset(map,chr<23) #Keep autosomal SNPs only (chr1-22)
# perform SNP matching
info_snp <- snp_match(sumstats, map)
# Assign the genotype to a variable for easier downstream analysis
genotype <- obj.bigSNP$genotypes

# Rename the data structures
CHR <- map$chr
POS <- map$pos

# get the CM information from 1000 Genome
# will download the 1000G file to the current directory (".")
POS2 <- snp_asGeneticPos(CHR, POS, dir = ".", ncores = NCORES)
# calculate LD
for (chr in 1:22) {
    # Extract SNPs that are included in the chromosome
    ind.chr <- which(info_snp$chr == chr)
    ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
    # Calculate the LD
    corr0 <- snp_cor(genotype, ind.col = ind.chr2, ncores = NCORES, infos.pos = POS2[ind.chr2], size = 3 / 1000)
    if (chr == 1) {
        ld <- Matrix::colSums(corr0^2)
        corr <- as_SFBM(corr0, tmp)
    } else {
        ld <- c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))
    }
}

# We assume the fam order is the same across different chromosomes
fam.order <- as.data.table(obj.bigSNP$fam)
# Rename fam order
setnames(fam.order,c("family.ID", "sample.ID"),c("FID", "IID"))


#Step 3: Perform LD score regression on betas
df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
ldsc <- snp_ldsc(ld,  length(ld), chi2 = (df_beta$beta / df_beta$beta_se)^2, sample_size = df_beta$n_eff, blocks = NULL)
h2_est <- ldsc[["h2"]]

#Create results DF for *.bimbam output
results<-data.frame()

#Step 4: Calculate null R^2 for binary trait
#library(rms)
# Reformat the phenotype file such that y is of the same order as the sample ordering in the genotype file
#y <- pheno[fam.order, on = c("FID", "IID")]
# use glm for binary trait (will also need the fmsb package to calculate the pseudo R2)
#null.formula <- paste("PCA", 1:6, sep = "", collapse = "+") %>% paste0("ADHD~Sex+Age+Carrier+", .) %>% as.formula 
#null.model<- lrm(null.formula, data = y)
#null.r2 <- null.model$stats["R2"]

#Step 5: Calculate PRS using infinitesemal model
#beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
# calculate PRS for all samples
#ind.test <- 1:nrow(genotype)
#pred_inf <- big_prodVec(genotype,beta_inf,ind.row = ind.test,ind.col = info_snp$`_NUM_ID_`)

#Test performance with regression model
#reg.formula <- paste("PCA", 1:6, sep = "", collapse = "+") %>% paste0("ADHD~PRS+Sex+Age+Carrier+", .) %>% as.formula
#reg.dat <- y
#reg.dat$PRS <- pred_inf
#inf.model <- lrm(reg.formula, dat=reg.dat)
#inf.r2 <- inf.model$stats["R2"]

#Scale PRS values (mean=1, SD=1) and reformat PRS to *.bimbam format for GEMMA
#prs_gemma<-scale(pred_inf)+1
#prs_gemma<-t(prs_gemma)
#prs_gemma<-cbind("ADHD_inf",0,0,prs_gemma)
#results<-rbind(results, prs_gemma)


##Step 6: Calculate PRS using grid model
## Prepare data for grid model
#p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2)
#h2_seq <- round(h2_est * c(0.7, 1, 1.4), 4)
#grid.param <-expand.grid(p = p_seq,h2 = h2_seq,sparse = c(FALSE, TRUE))
# Get adjusted beta from grid model
#beta_grid <-snp_ldpred2_grid(corr, df_beta, grid.param, ncores = NCORES)

# calculate PRS for all samples
#ind.test <- 1:nrow(genotype)
#pred_grid <- big_prodMat(genotype, beta_grid, ind.col = info_snp$`_NUM_ID_`)

#Loop through columns, scale PRS, and output in *.bimbam format
#for(i in 1:ncol(pred_grid)){
#    prs_grid <- pred_grid[,i]
#    prs_grid<-scale(pred_inf)+1
#    prs_grid<-t(prs_grid)
#    prs_grid<-cbind(paste0("ADHD_grid_",i),0,0,prs_grid)
#    results<-rbind(results, prs_grid)
#}

#Test performance with regression model
#reg.formula <- paste("PCA", 1:6, sep = "", collapse = "+") %>% paste0("ADHD~PRS+Sex+Age+Carrier+", .) %>% as.formula
#reg.dat <- y

#For below analyses, skip columns of the grid that fail with singularity error when generating LRM (not sure why this is happening)
# ~53% of columns fail this way
#grid.r2<-vector()
#index.pass<-vector()
#for(i in 1:ncol(pred_grid)){
#	flag<-TRUE
#    reg.dat$PRS <- pred_grid[,i]
#    grid.model <- tryCatch(lrm(reg.formula, dat=reg.dat),warning=function(w) flag<<-FALSE)
#    if (!flag) next
#    grid.r2 <- c(grid.r2, grid.model$stats["R2"])  
#    index.pass <- c(index.pass, i)
#    }

#max.r2 <- max(grid.r2)


#Step 7: Calculate PRS using auto model
# Get adjusted beta from the auto model (batch run needed for this step, run time ~5h on imputed data)
multi_auto <- snp_ldpred2_auto(corr,df_beta,h2_init = h2_est,vec_p_init = seq_log(1e-4, 0.9, length.out = NCORES),ncores = NCORES)
beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)

# calculate PRS for all samples
ind.test <- 1:nrow(genotype)
pred_auto <-big_prodMat(genotype,beta_auto,ind.row = ind.test,ind.col = info_snp$`_NUM_ID_`)
# scale the PRS generated from AUTO
pred_scaled <- apply(pred_auto, 2, sd)
final_beta_auto <-rowMeans(beta_auto[,abs(pred_scaled -median(pred_scaled)) < 3 * mad(pred_scaled)])
pred_auto <-big_prodVec(genotype,final_beta_auto,ind.row = ind.test,ind.col = info_snp$`_NUM_ID_`)

#Scale PRS values (mean=1, SD=1) and reformat PRS to *.bimbam format for GEMMA
prs_gemma<-scale(pred_auto)+1
prs_gemma<-t(prs_gemma)
prs_gemma<-cbind("educational_attainment_auto",0,0,prs_gemma)
results<-rbind(results, prs_gemma)

prs_auto<-t(pred_auto)
prs_auto<-cbind("educational_attainment_auto_noscale",0,0,prs_auto)
results<-rbind(results, prs_auto)

#Test performance with regression model
#reg.formula <- paste("PCA", 1:6, sep = "", collapse = "+") %>% paste0("ADHD~PRS+Sex+Age+Carrier+", .) %>% as.formula
#reg.dat <- y
#reg.dat$PRS <- pred_auto
#auto.model <- lrm(reg.formula, dat=reg.dat)
#auto.r2 <- auto.model$stats["R2"]

#Output results file as *.bimbam format
write.table(results,file="prs_gemma_educational_attainment.bimbam",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
