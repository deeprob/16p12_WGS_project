#!/bin/R




library(glue)
library(RareComb)



args = commandArgs(trailingOnly=T)
phenotype = args[1]
combos = as.numeric(args[2])

print(phenotype)
print(combos)

df = read.table('tables/variants.csv', sep=',', header=T)

# get list of cases and controls
cases = scan(glue('cases_controls/{phenotype}_cases.txt'), character(), quote = "")
controls = scan(glue('cases_controls/{phenotype}_controls.txt'), character(), quote = "")
all_samples = c(cases, controls)

# only keep samples that are cases or controls
df = df[df$Sample_Name %in% all_samples,]

# add column for output
df[,glue('Output_{phenotype}')] = as.integer(df$Sample_Name %in% cases)


combo_length=combos
max_freq_threshold=0.25
pval_filter_threshold=0.05
adj_pval_type='bonferroni'
min_power_threshold=0.7
sample_names_ind='Y'


if (combos == 2){
	min_indv_threshold = 5
} else {
	min_indv_threshold = 3
}

print(min_indv_threshold)


result = compare_enrichment(df,combo_length=combo_length, min_indv_threshold=min_indv_threshold, max_freq_threshold=max_freq_threshold, input_format='Input_', output_format='Output_', pval_filter_threshold=pval_filter_threshold, adj_pval_type=adj_pval_type, min_power_threshold=min_power_threshold, sample_names_ind=sample_names_ind)






filename = glue('statistics/3_rarecomb/{phenotype}_{combos}.csv')
write.csv(result, filename, row.names=F)






