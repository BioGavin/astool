#!/usr/bin/env python
# Usage: collect_gbk_file.py smiles.tsv target_folder
import os.path
import shutil
import sys
import pandas as pd


def get_file_path(tsv_file):
    info = pd.read_csv(tsv_file, sep='\t')
    for idx in info.index:
        folder_path = info.iloc[idx, 0]
        record_id = info.iloc[idx, 1]
        region_id = info.iloc[idx, 2]
        file_region_id = "000" + str(region_id)
        file_region_id = file_region_id[-3:]
        gbk_file_name = f'{record_id}.region{file_region_id}.gbk'
        gbk_file_path = os.path.join(folder_path, gbk_file_name)
        yield gbk_file_path


def mv2target(gbk, target_folder):
    changed_file_name = os.path.basename(os.path.dirname(gbk))[:15]
    target = os.path.join(target_folder, f'{changed_file_name}.{os.path.basename(gbk)}')
    shutil.copyfile(gbk, target)

if __name__ == '__main__':
    tsv_file, target_folder = sys.argv[1], sys.argv[2]
    gbk_files = get_file_path(tsv_file)
    for fp in gbk_files:
        mv2target(fp, target_folder)
