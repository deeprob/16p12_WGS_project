import sys
import csv

def genereturn(chrom,start,end,gene_lines):
	gene_list=[]
	for line in gene_lines:
		line=line.strip().split("\t")
		overlap=False
		gene_chrom=line[2]
		gene_start=int(line[4])
		gene_end=int(line[5])
		gene_name=line[12]
		if chrom==gene_chrom:
			if gene_start>=start and gene_start<=end: #Gene start within CNV
				overlap=True
			elif gene_end>=start and gene_end<=end: #Gene end within CNV
				overlap=True
			elif gene_start<=start and gene_end>=end: #CNV within gene
				overlap=True
		if overlap==True:
			gene_list.append(gene_name)
	return gene_list


#Open files
infile=open('output/merged.sv2.breakpoints_resolved.nearby_merged.tsv','r')
inlines=infile.readlines()
inlines = inlines[1:]

gene_file=open("/data4/software/annovar/humandb/hg19_refGene.txt","r")
gene_lines=gene_file.readlines()

outwriter = open("output/merged.sv2.genes.tsv",'w')
sample_id=""
for line1str in inlines:
	line=line1str.strip().split("\t")
	if line[0]!=sample_id:
		sample_id=line[0]
		print( "Processing sample "+sample_id)
	if line[2]=="chrom": #Header
		line.append("Genes")
		line.append("Gene_list")
	else:
		genes=genereturn(line[2],int(line[3]),int(line[4]),gene_lines)
		line.append(len(genes))
		if len(genes)>0:
			line.append(','.join(genes))
	line  = [str(x) for x in line]
	newline = '\t'.join(line)+'\n'
	outwriter.write(newline)
