#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=gangstr
#SBATCH -o logs/gangstr/submit_%a.log
#SBATCH -e logs/gangstr/submit_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=5G
#SBATCH --workdir /data5/16p12_WGS/structural_variants/vcf_callers
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user awt5304@psu.edu
#SBATCH --nodelist=qingyu
#SBATCH --array=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,73,75,76,77,79,80,82,83,84,87,88,89,90,91,95,96,98,99,100,101,102,105,106,107,108,109,110,111,112,116,117,118,122,123,127,128,129,130,132,133,134,135,136,137,138,140,141,142,143,145,148,149,150,151,153,154,156,157,160,161,162,163,164,165,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242

echo `date` starting job on $HOSTNAME as $USER

export SINGULARITY_CACHEDIR=/data5/16p12_WGS/structural_variants/vcf_callers/cache

fasta=/data/bx_references/hg19/ucsc.hg19.fasta
out_dir=/data5/16p12_WGS/structural_variants/vcf_callers/output/gangstr
tmp_dir=/data5/16p12_WGS/structural_variants/vcf_callers/tmp/gangstr

bamfile=`cat bam.list | awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id'`
sample=`echo $bamfile | cut -d"/" -f5`

echo `date` $bamfile
echo `date` $sample

singularity exec -B /data/bx_references/hg19,/data3/16p12_WGS,/data4/16p12_WGS,/data2/16p12_WGS \
	/data5/software/str-toolkit_latest.sif \
	GangSTR \
		--bam $bamfile \
		--ref $fasta \
		--regions hg19_ver13_1.bed \
		--out $tmp_dir/$sample \
		--include-ggl 

cp ${tmp_dir}/${sample}.vcf ${out_dir}/gangstr.${sample}.vcf
bgzip ${out_dir}/gangstr.${sample}.vcf
tabix -p vcf ${out_dir}/gangstr.${sample}.vcf.gz


echo `date` "done"





