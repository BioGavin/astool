#!/usr/bin/env python

import os
import sys

from Bio import SeqIO

from astool.antismash_utils import AntismashRegionGBKParser
from astool.utils import get_gbk_dir_ls

if __name__ == '__main__':
    input_gbk, outut_faa = sys.argv[1: 3]
    domain_record_ls = []
    for gbk_dir in get_gbk_dir_ls(input_gbk):
        gbk_file_name = os.path.basename(gbk_dir).rstrip('.gbk')
        gbk_file = AntismashRegionGBKParser(gbk_dir)
        for seq_rec in gbk_file.ex_a_domain():
            seq_rec.description += "fromfile: " + gbk_file_name
            domain_record_ls.append(seq_rec)
    SeqIO.write(domain_record_ls, outut_faa, "fasta")