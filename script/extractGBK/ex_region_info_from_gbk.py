#!/usr/bin/env python

import os.path
import sys

import pandas as pd

from astool.antismash_utils import AntismashRegionGBKParser


if __name__ == '__main__':
    gbk_dir_file, output_tsv = sys.argv[1], sys.argv[2]
    with open(gbk_dir_file, 'r') as f:
        gbk_dir_ls = f.read().splitlines()
    filedir, filename, bgc_type, bgc_length = [], [], [], []
    for gbk_dir in gbk_dir_ls:
        # print(gbk_dir)
        gbk_file = AntismashRegionGBKParser(gbk_dir)
        filedir.append(gbk_dir)
        filename.append(os.path.basename(gbk_dir))
        bgc_type.append(gbk_file.bgc_type)
        bgc_length.append(gbk_file.bgc_length)

    df = pd.DataFrame({
        'file_dir': filedir,
        'file_name': filename,
        'bgc_type': bgc_type,
        'bgc_length': bgc_length
    })
    df.to_csv(output_tsv, sep='\t', index=False)
