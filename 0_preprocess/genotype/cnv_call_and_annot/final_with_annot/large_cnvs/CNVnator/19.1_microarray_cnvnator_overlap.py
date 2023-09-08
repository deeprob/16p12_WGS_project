import pandas as pd
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Get counts of the overlaps
lookup_df = pd.read_csv('bed_files/19_overlap.bed', sep = '\t', header = None,
	names = ['Chr', 'Start', 'End', 'Type', 'Name', 'Length', 'Sample', 'Microarray_count', 'Intracohort_count',
		'microarray_freq', 'gnomADSV_AF', 'gnomADSV_AFfilter', 'gene_ids', 'gene_names', 'NEJM_Name', 'variant_id', 'inheritance',
		'micro_chr', 'micro_start', 'micro_end', 'micro_type', 'micro_sample'])

call_df = pd.read_csv('bed_files/17_inheritance.bed', sep = '\t')

# Function to find matching calls
def match(id):
	# Get all regions with 50% reciprocal overlap
	match_df = lookup_df[lookup_df.variant_id==id]
	# Get type matches
	match_df = match_df[((match_df.Type=='<DEL>') & (match_df.micro_type=='del')) | ((match_df.Type=='<DUP>') & (match_df.micro_type=='dup'))]
	if match_df.shape[0]==0:
		return False
	# Get sample matches
	match_df = match_df[match_df.Sample==match_df.micro_sample]
	if match_df.shape[0]==0:
		return False
	return True

call_df['match'] = call_df.variant_id.apply(match)
print(call_df.match.value_counts())

# Compare the QC metrics of the matching and non-matching calls

# Annotate each call with its QC metrics
# If a call is a merge of multiple calls, report all metrics
def add_metrics(row):
	sample = row['Sample']
	# Get call names
	names = row['Name'].split(';')

	# Get filename
	file = "/data5/16p12_WGS/structural_variants/sv_caller_postprocessing/2022_03_02/small_cnv_merge/compressed_caller_outputs/cnvnator/"+sample+".cnvnator.filtered.vcf.gz"

	# Get calls from the bgzipped VCF file
	lines = []
	for name in names:
		command = "bcftools view -i 'ID=\"%s\"' %s -H" % (name, file)
		match_lines = subprocess.run(command, capture_output = True, shell = True).stdout.decode()
		lines.append(match_lines)

	match_list = [i.split('\t') for i in lines if i!='']
	df = pd.DataFrame(match_list, columns = ['Chr','Pos', 'ID', 'Ref', 'Alt', 'Qual', 'FILTER', 'INFO', 'FORMAT', 'FORMAT_VALUES'])

	# We want to get:
	# 1. normalized_RD
	# 2. p-val1
	# 3. p-val2
	# 4. p-val3
	# 5. p-val4
	# 6. q0
	# Split these into their own columns
	df['norm_RD'] = df.INFO.apply(lambda info: float(info.split('natorRD=')[1].split(';')[0]))
	df['p1'] = df.INFO.apply(lambda info: float(info.split('natorP1=')[1].split(';')[0]))
	df['p2'] = df.INFO.apply(lambda info: float(info.split('natorP2=')[1].split(';')[0]))
	df['p3'] = df.INFO.apply(lambda info: float(info.split('natorP3=')[1].split(';')[0]))
	df['p4'] = df.INFO.apply(lambda info: float(info.split('natorP4=')[1].split(';')[0]))
	df['q0'] = df.INFO.apply(lambda info: float(info.split('natorQ0=')[1].split(';')[0]))

	# Return all values
	row['norm_RD'] = df['norm_RD'].to_list()
	row['p1'] = df['p1'].to_list()
	row['p2'] = df['p2'].to_list()
	row['p3'] = df['p3'].to_list()
	row['p4'] = df['p4'].to_list()
	row['q0'] = df['q0'].to_list()

	# "Best" value
	if np.mean(df['norm_RD']) > 1:
		row['best_norm_RD'] = max(df['norm_RD'])
	else:
		row['best_norm_RD'] = min(df['norm_RD'])
	row['best_p1'] = min(df['p1'])
	row['best_p2'] = min(df['p2'])
	row['best_p3'] = min(df['p3'])
	row['best_p4'] = min(df['p4'])
	row['best_q0'] = min(df['q0'])

	# "Worst" value
	if np.mean(df['norm_RD']) > 1:
		row['worst_norm_RD'] = min(df['norm_RD'])
	else:
		row['worst_norm_RD'] = max(df['norm_RD'])
	row['worst_p1'] = max(df['p1'])
	row['worst_p2'] = max(df['p2'])
	row['worst_p3'] = max(df['p3'])
	row['worst_p4'] = max(df['p4'])
	row['worst_q0'] = max(df['q0'])

	# Average value
	row['ave_norm_RD'] = np.mean(df['norm_RD'])
	row['ave_p1'] = np.mean(df['p1'])
	row['ave_p2'] = np.mean(df['p2'])
	row['ave_p3'] = np.mean(df['p3'])
	row['ave_p4'] = np.mean(df['p4'])
	row['ave_q0'] = np.mean(df['q0'])

	return row

anno_calls = call_df.apply(add_metrics, axis = 1)
print(anno_calls)

# Plot the values of "good" (matched) calls and "unknown" calls
# Add a new column for label
def labels(row):
	match = row['match']
	type = row['Type']

	if match and type=='<DEL>':
		return 'v_del'
	elif match and type=='<DUP>':
		return 'v_dup'
	elif not match and type=='<DEL>':
		return 'del'
	elif not match and type=='<DUP>':
		return 'dup'

anno_calls['Label'] = anno_calls.apply(labels, axis = 1)
print(anno_calls)
print(anno_calls.Label.value_counts())

# Make plots using best value, worst value, and average value
metrics = ['norm_RD', 'p1', 'p2', 'p3', 'p4', 'q0']
for super in ['best_', 'worst_', 'ave_']:
	fig, axs = plt.subplots(2, 3)
	plt.rc('xtick', labelsize=6)
	plt.rc('ytick', labelsize=6)
	row = 0
	col = 0
	for metric in metrics:
		axis = axs[row, col]
		sns.boxplot(x = 'Label', y = super+metric, data = anno_calls, ax = axis, fliersize = 0)
		sns.stripplot(x = 'Label', y = super+metric, data = anno_calls, ax = axis, color = 'k', size = 2)
		# Clean up the plots
		axis.set_title(metric)
		axis.set_ylabel(None)
		axis.set_xlabel(None)
		col+=1
		if col>2:
			col = 0
			row +=1

	plt.tight_layout()
	plt.savefig('qc_metrics_plots/'+super+'plots.pdf', dpi = 500)

# Just plot normalized read depth
fig, axs = plt.subplots(1, 3, sharey = True)
col = 0
for super in ['best_', 'worst_', 'ave_']:
	axis = axs[col]
	sns.boxplot(x = 'Label', y = super+'norm_RD', data = anno_calls, ax = axis, fliersize = 0)
	sns.stripplot(x = 'Label', y = super+'norm_RD', data = anno_calls, ax = axis, color = 'k', size = 2)
	axis.set_title(super.split('_')[0])
	axis.set_ylabel(None)
	axis.set_xlabel(None)
	axis.set(yscale='log')
	col+=1
plt.tight_layout()
plt.savefig('qc_metrics_plots/all_norm_RD.pdf', dpi = 500)
