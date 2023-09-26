#!/usr/bin/env python
# Usage: rename_gbk_in_metadata.py metadata new_metadata


import os
import sys

import pandas as pd


def change_gbk_name(metadata_dir):
    metadata = pd.read_csv(metadata_dir, sep='\t')
    metadata['gbk_name'] = None
    for idx in metadata.index:
        prefix = os.path.basename(metadata.iloc[idx, 0])[:15]
        record_id = metadata.iloc[idx, 1]
        region_id = metadata.iloc[idx, 2]
        file_region_id = "000" + str(region_id)
        file_region_id = file_region_id[-3:]
        old_gbk_file_name = f'{record_id}.region{file_region_id}'
        changed_gbk_file_name = f'{prefix}.{old_gbk_file_name}'
        metadata.iloc[idx, 4] = changed_gbk_file_name
    return metadata



if __name__ == '__main__':
    metadata_dir, new_metadata_dir = sys.argv[1], sys.argv[2]
    new_metadata = change_gbk_name(metadata_dir)
    new_metadata.to_csv(new_metadata_dir, sep='\t', index=False)