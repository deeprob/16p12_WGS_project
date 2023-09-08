import pandas as pd

# Combine large calls (CNVnator and microarray) and small calls into one file
cnvnator = pd.read_csv('/data5/16p12_WGS/structural_variants/sv_caller_postprocessing/2022_02_14/cnvnator/large_cnv_processing/final_calls/calls_by_gene_anno.txt', sep = '\t')
microarray = pd.read_csv('/data5/16p12_WGS/structural_variants/sv_caller_postprocessing/2022_03_09/annotate_penncnv/final_calls/calls_by_gene_anno.txt', sep = '\t')
small_cnv = pd.read_csv('/data5/16p12_WGS/structural_variants/sv_caller_postprocessing/2022_03_02/small_cnv_merge/final_calls/calls_by_gene_anno.txt', sep = '\t')

print(cnvnator)
print(microarray)
print(small_cnv)

df = pd.concat([microarray, small_cnv, cnvnator])
# Replace NaN with .
df.loc[df.NEJM.isnull(), 'NEJM']='.'
print(df)

# Remove duplicate genes
df2 = df.drop_duplicates(subset = ['Sample', 'Type', 'Gene_ID', 'Gene_Symbol'], keep = 'first')
print(df2)

# Look at duplicates, if you want
#print(df[df.duplicated(subset = ['Sample', 'Type', 'Gene_ID', 'Gene_Symbol'])])

df2.to_csv('sv_calls_combined.txt', sep = '\t', index = False)
