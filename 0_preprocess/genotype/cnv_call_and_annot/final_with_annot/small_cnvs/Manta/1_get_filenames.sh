#!/bin/bash

# Get filtered filenames
find /data5/16p12_WGS/structural_variants/vcf_callers/output/manta -name "*.manta.filtered.vcf" > manta_files.txt

# Add the SG code for easy filtering
python - << EOF

import pandas as pd

df = pd.read_csv('manta_files.txt', names = ['filename'])
df['Sample'] = df.filename.str.split('/', expand = True)[7].str.split('.', expand = True)[0]

df.to_csv('manta_files.csv', index = False)

EOF

rm manta_files.txt
