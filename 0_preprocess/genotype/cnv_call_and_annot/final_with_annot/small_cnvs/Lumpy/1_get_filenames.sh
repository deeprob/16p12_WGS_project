#!/bin/bash

# Get filtered filenames
find /data5/16p12_WGS/structural_variants/vcf_callers_2022_02_07/smoove_lumpy/vcfs/8_separate_samples/ -name "*_lumpy_filtered.vcf" > lumpy_files.txt

# Add the SG code for easy filtering
python - << EOF

import pandas as pd

df = pd.read_csv('lumpy_files.txt', names = ['filename'])
df['Sample'] = df.filename.str.split('/', expand = True)[8].str.split('_', expand = True)[0]

df.to_csv('lumpy_files.csv', index = False)

EOF

rm lumpy_files.txt
