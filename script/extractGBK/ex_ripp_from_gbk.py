#! /usr/bin/env python3
# Usage: ex_ripp_from_gbk.py gbk.list output.tsv

import os
import sys
import pandas as pd
from astool.antismash_utils import AntismashRegionGBKParser
from astool.utils import get_gbk_dir_ls

if __name__ == '__main__':
    gbk_input, tsv_output = sys.argv[1:3]
    gbk_name_ls = []
    leader_seqs = []
    core_seqs = []
    gbk_files = []
    get_gbk_dir_ls(gbk_input)
    for gbk_dir in get_gbk_dir_ls(gbk_input):
        gbk_file = AntismashRegionGBKParser(gbk_dir)
        leader_seq, core_seq = gbk_file.ex_lanthipeptide()
        for ls, cs in zip(leader_seq, core_seq):
            gbk_name_ls.append(os.path.basename(gbk_dir))
            leader_seqs.append(ls)
            core_seqs.append(cs)
            gbk_files.append(gbk_dir)

    pd.DataFrame({
        "gbk_file": gbk_files,
        "gbk_name": gbk_name_ls,
        "leader_seq": leader_seqs,
        "core_seq": core_seqs
    }).to_csv(tsv_output, sep='\t', index=False)
