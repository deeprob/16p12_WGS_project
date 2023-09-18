export PENNCNV_PATH=/data/software/PennCNV-1.0.5

#Input files needed: signal input text files for each sample with SNP name, B-allele freq., and log-R Ratio;
#.pfb (reference for allele B freq.), .hmm (expected numbers to build HMM--provided by software), and .gcmodel (%GC around each SNP) if available

#Step 1: Prepare input signal intensity files from raw Illumina Report
/data/software/PennCNV-1.0.5/split_illumina_report.pl -prefix raw_data/ -suffix .txt giriajan_grc_omni_2_190102_FinalReport.txt

#Step 2: Build .PFB file if not present from all signal input files. This will be located in the directory with the raw data files
#Download the snpposfile from Illumina (may need to manually remove "Heading" section and change assay heading from "MapInfo" to "Pos")
cut -d ',' -f 2,10,11 InfiniumOmniExpress-24v1-3_A1.csv > InfiniumOmniExpress-24v1-3_snp_list.csv
sed 's/,/\t/g' InfiniumOmniExpress-24v1-3_snp_list.csv > InfiniumOmniExpress-24v1-3_snp_list.tsv
$PENNCNV_PATH/compile_pfb.pl `ls raw_data/*.txt` --snpposfile InfiniumOmniExpress-24v1-3_snp_list.tsv --output omniexpress_2019.pfb

#Step 3A: Call CNVs per individual (autosomal)
perl $PENNCNV_PATH/detect_cnv.pl -test -hmm $PENNCNV_PATH/lib/hhall.hmm -pfb omniexpress_2019.pfb raw_data/*.txt -log penncnv_raw.log -out penncnv_raw.txt
#Step 3B: Call CNVs per individual on X chromosome. This uses a --sexfile argument that has the sex of each file (file name \t sex)
perl $PENNCNV_PATH/detect_cnv.pl -test --chrx -hmm $PENNCNV_PATH/lib/hhall.hmm -pfb omniexpress_2019.pfb raw_data/*.txt --sexfile omniexpress_2019_sexfile.tsv -log penncnv_raw_chrX.log -out penncnv_raw_chrX.tsv

#Note: check log file to see whether a "large fraction" of samples have WF value less than -0.04 or higher than 0.04--if this is the case, apply GCmodel adjustment
#Not applicable here (only 2 samples had |WF|>0.04)

#Step 4A: Call CNVs by family. Trio and quad tab-delimited files are in the format F M P (S). 
perl $PENNCNV_PATH/detect_cnv.pl -trio -hmm $PENNCNV_PATH/lib/hhall.hmm -pfb omniexpress_2019.pfb -cnv penncnv_raw.txt --listfile trios.txt -out penncnv_raw_family_trios.txt
perl $PENNCNV_PATH/detect_cnv.pl -quartet -hmm $PENNCNV_PATH/lib/hhall.hmm -pfb omniexpress_2019.pfb -cnv penncnv_raw.txt --listfile quads.txt -out penncnv_raw_family_quads.txt
#Step 4B: Call X-chromosome CNVs by family (split quad file into trios)
perl $PENNCNV_PATH/detect_cnv.pl -trio --chrx -hmm $PENNCNV_PATH/lib/hhall.hmm -pfb omniexpress_2019.pfb -cnv penncnv_raw_chrX.tsv --listfile trios_chrX.txt --sexfile omniexpress_2019_sexfile.tsv -out penncnv_raw_family_chrX_trios.txt

