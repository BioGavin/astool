import json
import os.path

import pandas as pd
import time
from Bio import SeqIO


def check_json(input_json) -> bool:
    """Check the validity of the input json file."""

    try:
        with open(input_json, 'r') as f:
            json.load(f)
        return True
    except ValueError:
        return False


def check_json_suffix(input_path) -> bool:
    """Check if the input file path has a json suffix."""
    if input_path.endswith("json"):
        return True
    else:
        return False


def check_gbk_suffix(input_path) -> bool:
    """Check if the input file path has a gbk suffix."""
    if input_path.endswith("gbk"):
        return True
    else:
        return False


def gen_dataframe(records):
    """Generate dataframe.

    Args:
        records: A generator containing named tuples.

    Returns: A Dataframe.

    """
    dataframe_ls = []
    for record_chunk in records:
        dataframe_ls.append(pd.DataFrame(record_chunk))
    return pd.concat(dataframe_ls)


def save_dataframe2tsv(df, fp):
    """Save dataframe to a tsv file.

    Args:
        df: dataframe.
        fp: tsv file path.

    Returns: None.

    """
    df.to_csv(fp, index=False, sep='\t')


def timer(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        print("%s running time: %s secs." % (func.__name__, end_time - start_time))
        return result

    return wrapper


def get_gbk_dir_ls(input_gbk):
    """Converting input to a list of paths as return."""
    gbk_dir_ls = []
    if check_gbk_suffix(input_gbk):
        gbk_dir_ls.append(input_gbk)
    else:
        with open(input_gbk, 'r') as f:
            gbk_dir_ls = f.read().splitlines()
    gbk_dir_ls = [gbk_dir.strip() for gbk_dir in gbk_dir_ls if check_gbk_suffix(gbk_dir.strip())]
    gbk_dir_ls = list(set(gbk_dir_ls))
    yield from gbk_dir_ls


def get_json_dir_ls(input_json):
    """Convert input to a list of paths as return."""
    json_dir_ls = []
    if check_json_suffix(input_json):
        json_dir_ls.append(input_json)
    else:
        with open(input_json, 'r') as f:
            json_dir_ls = f.read().splitlines()

    json_dir_ls = [json_dir.strip() for json_dir in json_dir_ls if check_json_suffix(json_dir.strip())]
    json_dir_ls = list(set(json_dir_ls))
    yield from json_dir_ls


if __name__ == '__main__':
    input_json = "/Users/zhouzhenyi/Documents/github/astool/test_data/json.txt"
    inp = "/Users/zhouzhenyi/Documents/github/astool/test_data/gbk.txt"
    gbk = "/Users/zhouzhenyi/Downloads/BGC0000001.gbk"
