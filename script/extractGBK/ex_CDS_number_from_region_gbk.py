import argparse

import pandas as pd

from astool.utils import get_gbk_dir_ls
from astool.antismash_utils import AntismashRegionGBKParser


def help():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
examples:
1. python ex_CDS_number_from_region_gbk.py -i gbk.list -o cds_number.tsv
        '''
    )

    parser.add_argument('-i', '--input', required=True, type=str,
                        help='Input a single region GBK file or a list of antiSMASH region GBK files.')
    parser.add_argument('-o', '--output', required=True, type=str,
                        help='Output TSV file.')

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = help()
    gbk_dirs = get_gbk_dir_ls(args.input)
    gbk_file_ls, cds_number_ls = [], []
    for gbk_file in gbk_dirs:
        gbk_file_ls.append(gbk_file)
        gbk = AntismashRegionGBKParser(gbk_file)
        cds_number = gbk.ex_cds_number()
        cds_number_ls.append(cds_number)
    output_df = pd.DataFrame({
        'gbk_file': gbk_file_ls,
        'cds_number': cds_number_ls
    })
    output_df.to_csv(args.output, sep='\t')
