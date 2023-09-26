#! /usr/bin/env python3

import os.path
import sys

import pandas as pd

from astool.antismash_utils import AntismashRegionGBKParser
from astool.utils import get_gbk_dir_ls

if __name__ == '__main__':
    gbk_input, tsv_output = sys.argv[1:3]

    gbk_name_ls = []
    region_completeness_ls = []

    for gbk_dir in get_gbk_dir_ls(gbk_input):
        gbk_name_ls.append(os.path.basename(gbk_dir))
        gbk_file = AntismashRegionGBKParser(gbk_dir)
        region_completeness = gbk_file.region_completeness
        region_completeness_ls.append(region_completeness)

    pd.DataFrame({
        "gbk_name": gbk_name_ls,
        "region_completeness": region_completeness_ls
    }).to_csv(tsv_output, sep='\t', index=False)
