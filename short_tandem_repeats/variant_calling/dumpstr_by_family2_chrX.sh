#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=dumpstr
#SBATCH -o logs/mergestr_by_family_chrX/dumpstr2.log
#SBATCH -e logs/mergestr_by_family_chrX/dumpstr2.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=2G
#SBATCH --chdir /data5/16p12_WGS/structural_variants/vcf_callers
#SBATCH --nodelist=durga

# from https://www.nature.com/articles/s41586-020-03078-7#Sec6:
# The merged VCF was then used as input to dumpSTR to compute locus-level filters using the parameters --min-locus-hwep 10−5 --min-locus-callrate 0.8 --filter-regions GRCh38GenomicSuperDup.sorted.gz --filter-regions-names SEGDUP to remove genotypes overlapping segmental duplications
#  For chromosome X, the Hardy–Weinberg equilibrium filter was computed based only on females. Filters obtained from analysing each phase were combined and any TRs failing locus-level filters in any phase were removed from further analysis.

echo `date` starting job on $HOSTNAME as $USER


export TMPDIR=/data5/SV2_tmp
export SINGULARITY_TMPDIR=/data5/SV2_tmp
export SINGULARITY_CACHEDIR=/data5/16p12_WGS/structural_variants/vcf_callers/durga_cache


vcf=output/mergestr_by_family_chrX/mergestr.vcf


singularity exec \
	-B /data5 \
	durga_cache/str-toolkit.simg \
	echo setting TMPDIR to /data5/sv2_temp; export TMPDIR=/data5/SV2_tmp; \
	dumpSTR \
	--vcf $vcf \
	--out output/mergestr_by_family_chrX/mergestr_filtered \
	--min-locus-callrate 0.8 \
	--filter-regions GRCh37genomicSuperDups.sorted.bed.gz  \
	--filter-regions-names SEGDUP


# Note: --min-locus-hwep 0.00001 is 10^-5
# Note: see get_GRCh37genomicSuperDups.sh to see where GRCh37genomicSuperDups.sorted.bed.gz comes from.

echo `date` finished

