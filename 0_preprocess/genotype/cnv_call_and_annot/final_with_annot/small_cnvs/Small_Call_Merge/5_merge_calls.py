import pandas as pd
import sys

# For each person, merge calls from all 3 callers into a combined callset
# Calls will only be merged on the basis of 50% reciprocal overlap
sample = sys.argv[1]
infile = '4_combined_calls/'+sample+'_combined.bed'
calls = pd.read_csv(infile, sep = '\t', names = ['Chr', 'Start', 'End', 'Type', 'ID', 'Length', 'Sample', 'Caller1_overlap', 'Caller2_overlap', 'Caller3_overlap'])

# The last 2 columns are not used anymore and can be removed
calls.drop(['Caller1_overlap', 'Caller2_overlap', 'Caller3_overlap'], axis = 1, inplace = True)
# The "Type" column needs to be unified, because different callers have different notations for Type
type_dict = {'<DEL>':'DEL', '<DUP>':'DUP', 'DEL':'DEL', 'DUP':'DUP'}
calls['Type'] = calls.Type.map(type_dict)
#print(calls)

# Get calls with a 50% reciprocal overlap
def recip_50(row):
	chr = row['Chr']
	start = row['Start']
	end = row['End']
	type = row['Type']
	length = row['Length']
	id = row['ID']

	# Find other calls that have a 50% reciprocal overlap with the current call
	same_chrom = calls[(calls.Chr==chr) & (calls.ID!=id) & (calls.Type==type)]
	overlap = same_chrom[((same_chrom.Start <= start) & (same_chrom.End >= start)) |
				((same_chrom.Start <= end) & (same_chrom.End >= end)) |
				((same_chrom.Start >= start) & (same_chrom.End <= end))]
	# If there are no other calls that meet these criteria
	# They may have been missed in earlier filtering because earlier filtering did not consider SV type
	if overlap.shape[0]==0:
		return False
	overlap['overlap50'] = overlap.apply(lambda s: float(max(end, s['End']) - min(start, s['Start']))/max(length, s['Length']) >= 0.5, axis = 1)
	overlap = overlap[overlap.overlap50]
	if overlap.shape[0]==0:
		return False
	# Return the IDs of the calls with a 50% reciprocal overlap
	ids = list(overlap.ID.unique())
	ids.sort()
	return ids

calls['direct_overlaps'] = calls.apply(recip_50, axis = 1)
calls = calls[calls.direct_overlaps!=False]
#print(calls)

# Combine list of 50% rciprocal overlap IDs
def combine_ids(row):
	in_ids = row['direct_overlaps']
	out_ids = in_ids

	# Iteratively add in extended overlaps until the list stops changing
	while True:
		for id in in_ids:
			add_ids = calls[calls.ID==id]['direct_overlaps']
			for a_id in add_ids:
				for b_id in a_id:
					if b_id not in out_ids:
						out_ids.append(b_id)

		if in_ids==out_ids:
			break

		in_ids = out_ids

	out_ids.sort()
	return(out_ids)

calls['indirect_overlaps'] = calls.apply(combine_ids, axis = 1)
#print(calls)

# Write the merge of the indirect overlap lists
def merge_indirect(indirect):
	to_merge = calls[calls.ID.isin(indirect)]

	# Get information for the final call
	chrom = ' '.join(list(to_merge.Chr.unique()))
	start = min(to_merge.Start)
	end = max(to_merge.End)
	type = ' '.join(list(to_merge.Type.unique()))
	id = ';'.join(indirect)
	length = end - start
	sample = ' '.join(list(to_merge.Sample.unique()))

	# Return merged call as series
	s = pd.Series([chrom, start, end, type, id, length, sample], index = ['Chr', 'Start', 'End', 'Type', 'ID', 'Length', 'Sample'])

	return s

# Only perform merge on unique sets of IDs
unique_sets = []
out_calls = []
for i in calls['indirect_overlaps']:
	if i not in unique_sets:
		unique_sets.append(i)
		out_calls.append(merge_indirect(i))

# Write final calls to file
merged_calls = pd.DataFrame(out_calls)
merged_calls.to_csv('merged_calls/'+sample+'_merged.bed', index = False, sep = '\t')
