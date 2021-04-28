import sys
import numpy as np

# load bed file

bedfile=open(sys.argv[1], 'r')
# bedfile=open('filter_regions.bed', 'r')

# format bed as numpy array
bed = bedfile.readlines()
bed = bed[1:]
bed = [s.strip().split('\t')[:3] for s in bed]
bed = np.array(bed)

beddir = {}
for i in list(range(1, 23)) + ['X', 'Y', 'M']:
	chrom = 'chr{}'.format(i)
	beddir[chrom] = bed[bed[:,0] == chrom][:,1:3].astype(int)


bedfile.close()


#Open SV vcf
infile=open(sys.argv[2], 'r')
# infile=open('/data5/SV2_tmp/filtered-cnv-calls/SG001.manta.filtered.dups_dels_only.vcf', 'r')
inlines=infile.readlines()

outfile=open(sys.argv[3], 'w')

for line in inlines:
	tabline=line.split('\t')
	if '#' in tabline[0]: #Gets all header lines
		outfile.write(line)
		pass
	else:
		chrom=tabline[0]
		start=int(tabline[1])
		info = tabline[7].split(';')
		
		# get end
		end=''
		for item in info:
			lst=item.split('=')
			if lst[0]=='END':
				end=int(lst[1])
		sv_length = end - start
		
		# check if overlaps
		bed = beddir[chrom]
		min_end   = np.minimum(bed[:,1], end)
		max_start = np.maximum(bed[:,0], start) 
		overlap = (min_end-max_start) > 0
		obed = bed[overlap]
		bp_that_overlap = (np.minimum(obed[:,1], end) - np.maximum(obed[:,0], start)).sum()
		
		is_overlapped = (bp_that_overlap/sv_length > 2/3)
				
		
		if is_overlapped:
			continue
		else:
			outfile.write(line)



infile.close()
outfile.close()
