#!/bin/bash

# Get filtered filenames
find /data5/16p12_WGS/structural_variants/vcf_callers/output/cnvnator/bin200 -name "cnv2vcf.*.cnvnator.filtered.vcf" > cnvnator_files.txt

# Add the SG code for easy filtering
python - << EOF

import pandas as pd

df = pd.read_csv('cnvnator_files.txt', names = ['filename'])
df['Sample'] = df.filename.str.split('/', expand = True)[8].str.split('.', expand = True)[1]

df.to_csv('cnvnator_files.csv', index = False)

EOF

rm cnvnator_files.txt
