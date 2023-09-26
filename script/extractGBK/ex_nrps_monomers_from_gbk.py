#!/usr/bin/env python

import os
import sys

import pandas as pd

from astool.antismash_utils import AntismashRegionGBKParser
from astool.utils import get_gbk_dir_ls

if __name__ == '__main__':
    input_gbk, output_tsv = sys.argv[1:3]
    file_name, locus_tag, monomers, count = [], [], [], []
    for gbk_dir in get_gbk_dir_ls(input_gbk):
        try:
            gbk_file = AntismashRegionGBKParser(gbk_dir)
            gbk_file_name = os.path.basename(gbk_dir).rstrip('.gbk')
            monomers_dict = gbk_file.get_nrps_pks_monomers()
            for k, v in monomers_dict.items():
                file_name.append(gbk_file_name)
                locus_tag.append(k)
                monomers.append('-'.join(v))
                count.append(len(v))
        except:
            print(gbk_dir)
    pd.DataFrame({
        "file_name": file_name,
        "locus_tag": locus_tag,
        "monomers": monomers,
        "count": count
    }).to_csv(output_tsv, sep='\t', index=False)