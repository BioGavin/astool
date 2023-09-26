#!/usr/bin/env python

import os.path
import sys

import pandas as pd

from astool.antismash_utils import AntismashRegionGBKParser


def ex_asdomain_number(gbk_dir):
    gbk_file = AntismashRegionGBKParser(gbk_dir)
    number_domain = len(gbk_file.ex_domain())
    return number_domain


if __name__ == '__main__':
    gbk_dir_file, output_tsv = sys.argv[1], sys.argv[2]
    with open(gbk_dir_file, 'r') as f:
        gbk_dir_ls = f.read().splitlines()
    filename, bgc_type, bgc_length, number_domain, smiles = [], [], [], [], []
    for gbk_dir in gbk_dir_ls:
        # print(gbk_dir)
        gbk_file = AntismashRegionGBKParser(gbk_dir)
        if not gbk_file.ex_smiles():
            smiles.append(None)
            filename.append(os.path.basename(gbk_dir))
            bgc_type.append(gbk_file.bgc_type)
            bgc_length.append(gbk_file.bgc_length)
            number_domain.append(len(gbk_file.ex_domain()))
        else:
            for smi in gbk_file.ex_smiles():
                smiles.append(smi)
                filename.append(os.path.basename(gbk_dir))
                bgc_type.append(gbk_file.bgc_type)
                bgc_length.append(gbk_file.bgc_length)
                number_domain.append(len(gbk_file.ex_domain()))
    df = pd.DataFrame({
        'filename': filename,
        'bgc_type': bgc_type,
        'bgc_length': bgc_length,
        'number_domain': number_domain,
        'smiles': smiles
    })
    df.to_csv(output_tsv, sep='\t', index=False)
