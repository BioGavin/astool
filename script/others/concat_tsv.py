#! /usr/bin/env python

import os
import sys

import pandas as pd


def concat_tsv(folder):
    files = os.listdir(folder)
    df_ls = []
    for f in files:
        df = pd.read_csv(os.path.join(folder, f), sep='\t')
        df_ls.append(df)
    concat_df = pd.concat(df_ls)
    return concat_df


def save2tsv(df, fn):
    df.to_csv(fn, index=False, sep='\t')


if __name__ == '__main__':
    folder, fn = sys.argv[1], sys.argv[2]
    concat_df = concat_tsv(folder)
    save2tsv(concat_df, fn)
