#!/bin/bash

# Get filtered filenames
find /data5/16p12_WGS/structural_variants/vcf_callers_2022_02_07/delly/vcfs/6_separate_samples/ -name "*_delly_filtered.vcf" > delly_files.txt

# Add the SG code for easy filtering
python - << EOF

import pandas as pd

df = pd.read_csv('delly_files.txt', names = ['filename'])
df['Sample'] = df.filename.str.split('/', expand = True)[8].str.split('_', expand = True)[0]

df.to_csv('delly_files.csv', index = False)

EOF

rm delly_files.txt
