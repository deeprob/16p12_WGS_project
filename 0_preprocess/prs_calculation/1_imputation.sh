#For imputation step--create two sets of data, one with all QC-passed data and one with European samples only
#Remove data for 7 families with non-European ancestry: PSU034 (microarray only), PSU055, PSU056, PSU061, PSU072, PSU080, and PSU083
/data5/software/plink --bfile ../16p12_microarray_dedup --make-bed --keep ../16p12_microarray_final_FAM_QC.valid.sample --out 16p12_microarray_impute_all --extract ../16p12_microarray_dedup.snplist
/data5/software/plink --bfile ../16p12_microarray_dedup --make-bed --keep 16p12_microarray_final_FAM_QC_EUR.valid.sample --out 16p12_microarray_impute_EUR --extract ../16p12_microarray_dedup.snplist

#Script adapted from https://imputationserver.readthedocs.io/en/latest/prepare-your-data/

#Create frequency file
/data5/software/plink --freq --bfile 16p12_microarray_impute_all --out 16p12_microarray_impute_all
/data5/software/plink --freq --bfile 16p12_microarray_impute_EUR --out 16p12_microarray_impute_EUR

#Run scripts to perform QC and split by chromosome
perl HRC-1000G-check-bim.pl -b 16p12_microarray_impute_all.bim -f 16p12_microarray_impute_all.frq -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h
sh Run-plink.sh

#Generate VCF files per chromosome
for i in {1..23}; do /data5/software/gotcloud/bin/vcfCooker vcfCooker --in-bfile 16p12_microarray_impute_all-updated-chr${i} --ref /data/bx_references/b37/human_g1k_v37.fasta --out 16p12_microarray_impute_all-updated-chr${i}.vcf --write-vcf; done
for i in {1..23}; do bgzip 16p12_microarray_impute_all-updated-chr${i}.vcf; tabix -p vcf 16p12_microarray_impute_all-updated-chr${i}.vcf.gz; done
for i in {1..23}; do bcftools sort 16p12_microarray_impute_all-updated-chr${i}.vcf.gz -Oz -o 16p12_microarray_impute_all-updated-sort-chr${i}.vcf.gz; done
#Optional QC check
for i in {1..23}; do python /data5/software/checkVCF.py -r /data/bx_references/b37/human_g1k_v37.fasta -o out 16p12_microarray_impute_all-updated-sort-chr${i}.vcf.gz; done

#After running imputation--merge VCFs by chromosome
bcftools concat -o 16p12_microarray_TOPMED_impute.vcf.gz -Oz *dose.vcf.gz

#Create BED/BIM/FAM file from imputed VCF
/data5/software/plink --vcf 16p12_microarray_TOPMED_impute.vcf.gz --out 16p12_microarray_TOPMED_impute

#FAM file is identical to before, so OK to copy previous file
rm 16p12_microarray_TOPMED_impute.fam
cp ../16p12_microarray_impute_all.fam 16p12_microarray_TOPMED_impute.fam

#Rerun QC steps from initial pipeline

#PLINK QC filters: Cohort MAF <0.05, Hardy-Weinberg equilibrium <1e-6 (under selection), geno<0.01 (SNPs missing in >1% of subjects), mind<0.01 (Individuals with >1% missing genotypes)
/data5/software/plink --bfile 16p12_microarray_TOPMED_impute --maf 0.05 --hwe 1e-6 --geno 0.01 --make-bed --out 16p12_microarray_final_QC
/data5/software/plink --bfile 16p12_microarray_final_QC --mind 0.01 --write-snplist --make-bed --out 16p12_microarray_final_FAM_QC
#5,691,763 variants retained; no individual samples removed (8 samples were removed from previous pipeline)

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
#5 samples were removed: SG173, SG201, SG230, SG268, SG402
#Also remove poor quality and 15q13.3 samples--353 samples remain

#Check sex of samples: Generate sexcheck file in PLINK
/data5/software/plink --bfile 16p12_microarray_final_FAM_QC --extract 16p12_microarray_final_FAM_QC.prune.in --keep 16p12_microarray_final_FAM_QC.valid.sample --check-sex --out 16p12_microarray_final_FAM_QC
#All samples have matching sexes--so skip the removal step

#Calculate principal component eigenvectors for population stratification
/data5/software/plink --bfile 16p12_microarray_final_FAM_QC --extract 16p12_microarray_final_FAM_QC.prune.in --pca 6 --out 16p12_microarray_final_FAM_QC

#Remove duplicate SNPs from PLINK files
cut -f 2 16p12_microarray_final_FAM_QC.bim | sort | uniq -d > 16p12_dup_snps.txt
#No duplicate SNPs identified, so skipped step

#Output final general file. Exclude SNPs within hg38 16p12.1 deletion/adjacent region (chr16:21788679-22588679); 5,497,394 SNPs in final set
/data5/software/plink --bfile 16p12_microarray_final_FAM_QC --make-bed --keep 16p12_microarray_final_FAM_QC.valid.sample --out 16p12_microarray_final --extract 16p12_microarray_final_FAM_QC.snplist --write-snplist --exclude-snps chr16:21783909:G:A-chr16:22701658:G:T

#Lift over PLINK file to hg19 from hg38 (note: done at end to reduce runtime post-filtering; convert to and from .ped/.map format)
/data5/software/plink --bfile 16p12_microarray_final --recode tab --out 16p12_microarray_final
#Run liftover script (https://github.com/sritchie73/liftOverPlink/blob/master/liftOverPlink.py; **use Python2**)
python2 liftOverPlink.py -m 16p12_microarray_final.map -o 16p12_microarray_final_hg19 -p 16p12_microarray_final.ped -c hg38ToHg19.over.chain.gz -e /data/software/liftOver
/data5/software/plink --file 16p12_microarray_final_hg19 --allow-extra-chr --chr 1-22 XY --make-bed --write-snplist --recode --out 16p12_microarray_final_hg19

#Copy BIM file to retain original for downstream strand-flipping analysis
cp 16p12_microarray_final_hg19.bim 16p12_microarray_final_hg19_original.bim

#Adjust mismatched SNPs using strand-flipping in R
#See file PRS_file_processing.R


