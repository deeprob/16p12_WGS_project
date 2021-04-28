#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=mantassd
#SBATCH -o logs/manta/submit_%a.log
#SBATCH -e logs/manta/submit_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=2G
#SBATCH --workdir /data5/16p12_WGS/structural_variants/vcf_callers
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user awt5304@psu.edu
#SBATCH --nodelist=qingyu
#SBATCH --array=1-205%1

echo `date` starting job on $HOSTNAME as $USER

fasta=/data/bx_references/hg19/ucsc.hg19.fasta
out_dir=/data5/16p12_WGS/structural_variants/vcf_callers/output/manta
tmp_dir=/data5/16p12_WGS/structural_variants/vcf_callers/tmp/manta
bedfile=/data5/16p12_WGS/structural_variants/vcf_callers/chroms.bed.gz

bamfile=`cat bam.list | awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id'`
sample=`echo $bamfile | cut -d"/" -f5`

echo `date` $bamfile
echo `date` $sample

manta_dir=/data3/software/manta-1.6.0.centos6_x86_64

mkdir -p ${tmp_dir}/$sample

${manta_dir}/bin/configManta.py \
	--bam=${bamfile} \
	--referenceFasta=${fasta} \
	--runDir=${tmp_dir}/${sample} \
	--callRegions=$bedfile

# theres a bug with the line isEmail = isLocalSmtp()
# so replace that with 		 isEmail = False
cat ${tmp_dir}/${sample}/runWorkflow.py | sed 's/isEmail = isLocalSmtp()/isEmail = False/g' > ${tmp_dir}/${sample}/runWorkflow.noemail.py
chmod +x ${tmp_dir}/${sample}/runWorkflow.noemail.py
${tmp_dir}/${sample}/runWorkflow.noemail.py -j 40 -g 40

gunzip --keep ${tmp_dir}/${sample}/results/variants/diploidSV.vcf.gz
mv ${tmp_dir}/${sample}/results/variants/diploidSV.vcf ${out_dir}/${sample}.manta.vcf

echo `date` "done"




