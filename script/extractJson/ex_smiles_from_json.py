#!/usr/bin/env python
# Usage: ex_smiles json_dir output.tsv

import sys


from astool.ex_smiles import get_antismash_smiles_records
from astool.utils import gen_dataframe, save_dataframe2tsv


if __name__ == '__main__':
    # json_dir, output, type, verbose = args.json_dir, args.output, args.type, args.verbose
    json_dir, output = sys.argv[1: 3]
    print(json_dir, output)
    smiles_records = get_antismash_smiles_records(json_dir, verbose=True)
    smiles_df = gen_dataframe(smiles_records)
    save_dataframe2tsv(smiles_df, output)