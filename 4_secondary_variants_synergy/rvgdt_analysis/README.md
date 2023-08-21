# Describes all the scripts and their original locations required to run rvgdt

## Scripts available in dropbox
1. 1_prepare_sample_list.py

    Original Location: Dropbox/16p12.2 project/Human patients project/WGS paper/25_TDT tests/RV-GDT/1_prepare_sample_list.py

2. 2_keep_only_complete_trios.py

    Original Location: Dropbox/16p12.2 project/Human patients project/WGS paper/25_TDT tests/RV-GDT/3_keep_only_complete_trios.py

3. 3_remove_unkowns_from_ped.py

    Original Location: Dropbox/16p12.2 project/Human patients project/WGS paper/25_TDT tests/RV-GDT/4_remove_unkowns_from_ped.py

4. 4a_preprocess_for_rvgdt.sh

    Original Location: Dropbox/16p12.2 project/Human patients project/WGS paper/25_TDT tests/RV-GDT/5_rvgdt.sh

5. 4b_run_rvgdt.sh

    Original Location: Dropbox/16p12.2 project/Human patients project/WGS paper/25_TDT tests/RV-GDT/5_rvgdt.sh

6. 5_statistics.R

    Original Location: Dropbox/16p12.2 project/Human patients project/WGS paper/25_TDT tests/RV-GDT/7_statistics.R

## Scripts not available in dropbox
1. utils/sv2_split.py

    Original Location: /data5/anastasia/rvgdt/scripts/sv2_split.py

2. utils/gt_to_rvgdt.py

    Original Location: /data5/anastasia/rvgdt/scripts/gt_to_rvgdt.py

3. utils/map_phenotype_split_2.py

    Original Location: /data5/anastasia/rvgdt/scripts/map_phenotype_split_2.py


## Software
1. RVGDT
    Original Location: https://github.com/statgenetics/rv-gdt

## Environments
Note: Environment creation details are present in project root README.md file

1. Scripts that run under "rvgdt_pro_py" conda environment
    - 1_prepare_sample_list.py
    - 2_keep_only_complete_trios.py
    - 3_remove_unkowns_from_ped.py
    - 4a_preprocess_for_rvgdt.sh

2. Scripts that run under "rvgdt_python" conda environment
    - 4b_run_rvgdt.sh

3. Scripts that run under "rvgdt_pro_r" conda environment
    - 5_statistics.R
