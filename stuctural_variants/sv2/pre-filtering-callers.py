import sys

#Open SV vcf
infile=open(sys.argv[1], 'r')
inlines=infile.readlines()

outfile=open(sys.argv[2], 'w')

for line in inlines:
	tabline=line.split('\t')
	if '#' in tabline[0]: #Gets all header lines
		outfile.write(line)
	else:
		start=int(tabline[1])
		info = tabline[7].split(';')
		end=''
		for item in info:
			lst=item.split('=')
			if lst[0]=='END':
				end=int(lst[1])
		if end=='' or end-start<=10000000: #Removes very long SVs - these slow down SV2 and are likely artifacts anyway
			outfile.write(line)

infile.close()
outfile.close()
