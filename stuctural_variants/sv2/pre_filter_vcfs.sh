#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=filter
#SBATCH -o logs/pre_filter_%a.log
#SBATCH -e logs/pre_filter_%a.err
#SBATCH --ntasks=1
#SBATCH --time=720:00:00
#SBATCH --mem-per-cpu=50G
#SBATCH --chdir /data5/16p12_WGS/structural_variants/sv2
#SBATCH --nodelist=qingyu
#SBATCH --array=1-344

echo `date` start

#Change the tmp directory - using default directory causes issues
export TMPDIR=/data5/SV2_tmp

#Get SG code
code=$( head -$SLURM_ARRAY_TASK_ID all_codes.list | tail -1 )

echo `date` $code

# filter out calls that are greater than 10 Mbp
python3 pre-filtering-callers.py /data5/16p12_WGS/structural_variants/vcf_callers/output/manta/$code.manta.vcf vcfs/filtered/$code.manta.long_svs_filtered.vcf
python3 pre-filtering-callers.py /data5/16p12_WGS/structural_variants/vcf_callers/output/delly/$code.delly.vcf vcfs/filtered/$code.delly.long_svs_filtered.vcf
python3 pre-filtering-callers.py /data3/16p12_WGS/parsing_cnv_callers/sv2/lumpy_VCFs/$code.lumpy.vcf vcfs/filtered/$code.lumpy.long_svs_filtered.vcf
python3 pre-filtering-callers.py /data5/16p12_WGS/structural_variants/vcf_callers/output/cnvnator/bin200/cnv2vcf.$code.cnvnator.vcf vcfs/filtered/$code.cnvnator.long_svs_filtered.vcf

# filter calls that aren't for dups or dels
python3 filter_dups_dels.py vcfs/filtered/$code.manta.long_svs_filtered.vcf vcfs/filtered/$code.manta.long_svs_filtered.dups_dels_only.vcf
python3 filter_dups_dels.py vcfs/filtered/$code.delly.long_svs_filtered.vcf vcfs/filtered/$code.delly.long_svs_filtered.dups_dels_only.vcf
python3 filter_dups_dels.py vcfs/filtered/$code.lumpy.long_svs_filtered.vcf vcfs/filtered/$code.lumpy.long_svs_filtered.dups_dels_only.vcf
python3 filter_dups_dels.py vcfs/filtered/$code.cnvnator.long_svs_filtered.vcf vcfs/filtered/$code.cnvnator.long_svs_filtered.dups_dels_only.vcf


# filter calls with > 2/3 overlap with regions that are difficult to map
EXCLUDE_BED=filter_regions.bed

python3 filter_segdups.py $EXCLUDE_BED vcfs/filtered/$code.manta.long_svs_filtered.dups_dels_only.vcf vcfs/filtered/$code.manta.long_svs_filtered.dups_dels_only.segdups_filtered.vcf
python3 filter_segdups.py $EXCLUDE_BED vcfs/filtered/$code.delly.long_svs_filtered.dups_dels_only.vcf vcfs/filtered/$code.delly.long_svs_filtered.dups_dels_only.segdups_filtered.vcf
python3 filter_segdups.py $EXCLUDE_BED vcfs/filtered/$code.lumpy.long_svs_filtered.dups_dels_only.vcf vcfs/filtered/$code.lumpy.long_svs_filtered.dups_dels_only.segdups_filtered.vcf
python3 filter_segdups.py $EXCLUDE_BED vcfs/filtered/$code.cnvnator.long_svs_filtered.dups_dels_only.vcf vcfs/filtered/$code.cnvnator.long_svs_filtered.dups_dels_only.segdups_filtered.vcf


echo `date` end



