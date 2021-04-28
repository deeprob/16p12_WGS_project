#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=MPC
#SBATCH -o logs/MPC_%a.log
#SBATCH -e logs/MPC_%a.err
#SBATCH --cpus-per-task=2
#SBATCH --time=240:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --chdir /data5/SSC_WGS/SFARI_SSC_WGS_2/annotations
#SBATCH --array 1-24

echo `date` started on $HOSTNAME


# annotations for chroms ending with $a
a=$SLURM_ARRAY_TASK_ID

# if s > 22 set x,y,m chromosomes
if [[ $a -eq 23 ]]
then
  a=X
elif [[ $a -eq 24 ]]
then
  a=Y
fi

echo `date` working on chr${a}
echo `date` getting protein coding annotations

# get protein coding annotations
# python3 get_protein_coding.py vcfs/genes/chr${a}.hg38_multianno.vcf.gz | bgzip > vcfs/protein_coding/chr${a}.vcf.gz


echo `date` getting MPC

/data5/software/slivar  expr --vcf vcfs/protein_coding/chr${a}.vcf.gz \
	-g /data5/anastasia/liftover/mpc_hg38_gnotate.zip \
	| bgzip \
	> vcfs/protein_coding/chr${a}.with_MPC.vcf.gz




echo `date` getting sample genotypes

# to table and add samps
vcf_with_annotations=vcfs/protein_coding/chr${a}.with_MPC.vcf.gz
vcf_with_samples=vcfs/norm/chr${a}.vcf.gz
error_message_filename=vcfs/protein_coding/chr${a}_missed_alleles.txt
outfile=vcfs/protein_coding/chr${a}.txt
python3 add_sample_info.py $vcf_with_annotations $vcf_with_samples $error_message_filename > $outfile


echo `date` fin



