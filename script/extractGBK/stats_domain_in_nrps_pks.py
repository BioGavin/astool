#! /usr/bin/env python3

import os.path
import sys

import pandas as pd

from astool.antismash_utils import AntismashRegionGBKParser
from astool.utils import get_gbk_dir_ls


def get_nrps_pks_domain_name_df(gbk_dir):
    gbk_name = os.path.basename(gbk_dir)
    gbk_file = AntismashRegionGBKParser(gbk_dir)
    domain_name_ls = []
    for d in gbk_file.get_nrps_pks_domains():
        domain_name = "+".join(d.qualifiers.get('aSDomain'))
        domain_name_ls.append(domain_name)
    domain_name_df = pd.DataFrame({
        "gbk_name": gbk_name,
        "domain_name": domain_name_ls
    })
    return domain_name_df


if __name__ == '__main__':
    input_gbk, output_tsv = sys.argv[1:3]
    gbk_dir_ls = get_gbk_dir_ls(input_gbk)
    domain_name_df_ls = []
    for gbk_dir in gbk_dir_ls:
        domain_name_df_ls.append(get_nrps_pks_domain_name_df(gbk_dir))
    domain_name_df = pd.concat(domain_name_df_ls)
    res_df = domain_name_df.value_counts().unstack()
    res_df.to_csv(output_tsv, sep='\t')
    # print(domain_name_df)

    # test
    # gbk_dir = "/Users/zhouzhenyi/Documents/github/astool/test_data/antiSMASH/GCF_000196475.1_ASM19647v1_genomic/NC_012962.1.region016.gbk"
    # domain_name_df = get_nrps_pks_domain_name_df(gbk_dir)
    # print(domain_name_df)
    # print(domain_name_df.value_counts().unstack())
