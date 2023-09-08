
# Delly reports multiple SV types
# We only want deletions and duplications

import pandas as pd
import sys

file = sys.argv[1]
sample = sys.argv[2]

# Get index of header lines to skip
exclude = [i for i, line in enumerate(open(file, 'r')) if line.startswith('#')]

# Read file as dataframe
vcf = pd.read_csv(sys.argv[1], sep = '\t', skiprows = exclude, names = ['Chr','Pos', 'ID', 'Ref', 'Alt', 'Qual', 'FILTER', 'INFO', 'FORMAT', 'FORMAT_VALUES'])

# Make a new column for Type for easy filtering
vcf['Type'] = vcf.INFO.str.split('SVTYPE=', expand = True)[1].str.split(';', expand = True)[0]

# Remove calls that are not deletions or duplications
vcf[(vcf.Type=='DEL') | (vcf.Type=='DUP')].to_csv('vcfs/2_get_cnvs/'+sample+'_dels_dups.vcf', sep = '\t', index = False)
