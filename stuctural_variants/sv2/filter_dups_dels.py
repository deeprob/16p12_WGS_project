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
		svtype=''
		for item in info:
			lst=item.split('=')
			if lst[0]=='SVTYPE':
				svtype=(lst[1])
		if svtype == 'DUP' or svtype == 'DEL': #Removes anything that's not a dup or del
			outfile.write(line)

infile.close()
outfile.close()
