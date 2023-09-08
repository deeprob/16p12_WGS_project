# Get de novo calls
head -1 call_tables/6_genic_filter.txt > call_tables/7_denovo_calls.txt
grep 'de novo' call_tables/6_genic_filter.txt >> call_tables/7_denovo_calls.txt

# Visually confirm the de novo calls using the visualizations Matt had made
# Visualizations were copied from:
	# Dropbox\16p12.2 project\Human patients project\WGS paper\6_Microarray calls\Visualization\Second_hits_all_batches, OR
	# Dropbox\16p12.2 project\Human patients project\Exome and chips data\SNP arrays\Final_calls\Plots\probands
# Visualizations for putative de novo calls were compiled on Dropbox: Dropbox\16p12.2 project\Human patients project\WGS paper\6_Microarray calls\Inheritance\Denovo_plots
