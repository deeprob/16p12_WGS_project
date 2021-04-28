# running SV2

Scripts to run SV2

### pre_filter_vcfs.sh

Runs 3 python scripts to filter problematic calls from manta, lumpy, delly, and cnvnator.

1. pre-filtering-callers.py filters calls > 10 Mbp
2. filter_dups_dels.py filters calls that aren't dups or dels
3. filter_segdups.py filters segdups that are in the filter_regions.bed file

### sv2_submit.sh

Calls SV2 on filtered calls from manta, lumpy, delly, and cnvnator.

Note: for samples SG043, SG044, and SG138, I replaced the line 

    SNV_VCF="/data3/16p12_WGS/parsing_cnv_callers/sv2/WGS_VCFs/"$code"_SNVs.vcf.gz"

with 

    SNV_VCF=/data4/16p12_WGS/batch4_1019/SG138/SG138_gvcf.g.vcf.gz
    SNV_VCF=/data3/16p12_WGS/batch3_0618/SG043/SG043_gvcf.g.vcf.gz
    SNV_VCF=/data3/16p12_WGS/batch3_0618/SG044/SG044_gvcf.g.vcf.gz

Because of the error:

    FATAL ERROR: SNV file(s) were not formatted correctly

Second note: Sample SG231 wasn't in the all_codes.list file, so I ran that seperately.

