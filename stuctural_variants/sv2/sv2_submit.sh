#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=sv2
#SBATCH -o logs/sv2_%a.log
#SBATCH -e logs/sv2_%a.err
#SBATCH --ntasks=1
#SBATCH --time=720:00:00
#SBATCH --mem-per-cpu=50G
#SBATCH --chdir /data5/16p12_WGS/structural_variants/sv2
#SBATCH --nodelist=qingyu
#SBATCH --array=52,53,54,55,56,57,58,59,60,61,63,68,69,71,72,73,75,76,77,87,88,90,94,95,97,98,99,100,101,102,103,104,105,106,107,109,110,111,113,185,186,187,188,189,190,191,192,193,194,195,196,197,198,201,202,203,204,205,206,207,208,209,210,211,232,235,284,285,286,291,300,308,310,313,316,318,320,323,328,329,330,331,332,333,334,335,336,337

echo `date` start

#Change the tmp directory - using default directory causes issues
export TMPDIR=/data5/SV2_tmp

#Get SG code
code=$( head -$SLURM_ARRAY_TASK_ID all_codes.list | tail -1 )

echo `date` $code

PED=/data3/16p12_WGS/parsing_cnv_callers/sv2/all_batches.ped

BAM=$( grep $code /data3/16p12_WGS/parsing_cnv_callers/sv2/WGS_bam_list.list )
MANTA=/data5/16p12_WGS/structural_variants/sv2/vcfs/filtered/$code.manta.long_svs_filtered.dups_dels_only.segdups_filtered.vcf
DELLY=/data5/16p12_WGS/structural_variants/sv2/vcfs/filtered/$code.delly.long_svs_filtered.dups_dels_only.segdups_filtered.vcf
LUMPY=/data5/16p12_WGS/structural_variants/sv2/vcfs/filtered/$code.lumpy.long_svs_filtered.dups_dels_only.segdups_filtered.vcf
CNVNATOR=/data5/16p12_WGS/structural_variants/sv2/vcfs/filtered/$code.cnvnator.long_svs_filtered.dups_dels_only.segdups_filtered.vcf
SNV_VCF="/data3/16p12_WGS/parsing_cnv_callers/sv2/WGS_VCFs/"$code"_SNVs.vcf.gz"


PREP_FILE=/data5/16p12_WGS/structural_variants/sv2/sv2_preprocessing/
FEATS_FILE=/data5/16p12_WGS/structural_variants/sv2/sv2_features/
TMP_FILE=/data5/16p12_WGS/structural_variants/sv2/sv2_tmpfiles/$code

mkdir $TMP_FILE

sv2 -i $BAM -v $MANTA $LUMPY $DELLY $CNVNATOR -snv $SNV_VCF -p $PED -L sv2_logs/${code}.log -T $TMP_FILE -o ${code} -merge -min-ovr 0.5 # -pre $PREP_FILE

echo `date` end



