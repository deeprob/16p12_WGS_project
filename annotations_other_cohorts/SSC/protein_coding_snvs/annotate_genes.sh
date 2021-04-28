#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=genes
#SBATCH -o logs/annovar_%a.log
#SBATCH -e logs/annovar_%a.err
#SBATCH --cpus-per-task=2
#SBATCH --time=240:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --chdir /data5/SSC_WGS/SFARI_SSC_WGS_2/annotations
#SBATCH --array 1-24

echo `date` started on $HOSTNAME

export PATH=$PATH:/data5/software

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

vcf=vcfs/gnomad/chr${a}.vcf.gz
out=vcfs/genes/chr${a}

echo `date` working on chr${a}
echo `date` invcf=$vcf
echo `date` outvcf=$out

perl /data4/software/annovar/table_annovar.pl $vcf /data4/software/annovar/humandb/ \
	 -buildver hg38 \
	 -out $out \
	 -remove \
	 -protocol refGene \
	 -operation gx \
	 -nastring . \
	 -thread 6 \
	 -vcfinput \
	 -xref /data4/software/annovar/humandb/gene_annotations.txt \
	 -arg '-hgvs'


bgzip ${out}.hg38_multianno.vcf
tabix -p vcf ${out}.hg38_multianno.vcf

# cleanup
rm ${out}.avinput

echo `date` fin


