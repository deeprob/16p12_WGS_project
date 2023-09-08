#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=smoove
#SBATCH -o logs/7_filter.log
#SBATCH -e logs/7_filter.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=15G
#SBATCH --chdir /data5/16p12_WGS/structural_variants/vcf_callers_2022_02_07/smoove_lumpy
#SBATCH --nodelist qingyu


echo `date` starting job on $HOSTNAME as $USER


in_file=/data5/16p12_WGS/structural_variants/vcf_callers_2022_02_07/smoove_lumpy/vcfs/4_paste/merged.smoove.square.vcf.gz 
out_dir=/data5/16p12_WGS/structural_variants/vcf_callers_2022_02_07/smoove_lumpy/vcfs/7_filter/



export PATH=/data5/anastasia/sw/bcftools-1.12:$PATH


# Instruction at https://github.com/brentp/smoove
# There is also more inoformation about which filters to use at 
# https://github.com/brentp/smoove/issues/109

# brentp uses the following filters

# some notes on his command:
# LQHET stands for low quality heterozygous variant
# HQHET stands for high quality heterozygous variant
# HQHA stands for high quality homozygous variants

# slivar expr \
#   --info "variant.call_rate > 0.5 && ((INFO.SVTYPE == 'DEL') ||
# (INFO.SVTYPE == 'DUP'))" \
#   --sample-expr \
#     "LQHET:sample.alts == 1 && (((sample.DHFFC > 0.75) && (INFO.SVTYPE
# == 'DEL')) || ((sample.DHFFC < 1.25) && INFO.SVTYPE == 'DUP'))" \
#   --sample-expr \
#     "HQHET:sample.alts == 1 && (((sample.DHFFC < 0.75) && (INFO.SVTYPE
# == 'DEL')) || ((sample.DHFFC > 1.25) && INFO.SVTYPE == 'DUP'))" \
#   --sample-expr \
#     "HQHA:sample.alts == 2 && (((sample.DHFFC < 0.5) && (INFO.SVTYPE
# == 'DEL')) || ((sample.DHFFC > 1.5) && INFO.SVTYPE == 'DUP'))" \
#   -o $out \
#   -p $ped \
#   --vcf $vcf

# Only keep heterozygous dels with DHFFC < 0.75, heterozygous dups with DHFFC > 1.25, homozygous dels with DHFFC < 0.5, and homozygous dups with DHFFC > 1.5
bcftools filter --set-GT . -i '(DHFFC < 0.75 & SVTYPE="DEL" & GT="HET" & GT="ALT") |
	(DHFFC > 1.25 & SVTYPE="DUP" & GT="HET" & GT="ALT") |
	(DHFFC < 0.5 & SVTYPE="DEL" & GT="HOM" & GT="ALT") |
	(DHFFC > 1.5 & SVTYPE="DUP" & GT="HOM" & GT="ALT")' $in_file > $out_dir/merged.filtered.vcf



# view possible GT fields
bcftools query -f '[%GT\t%DHFFC\t%SVTYPE\n]' $out_dir/merged.filtered.vcf | cut -f1 | sort | uniq -c






echo `date` finished




