#!/usr/bin/env python
# Usage: download_antismash_db.py tsv target_path


import os.path
import re
import sys
from multiprocessing import Pool
import pandas as pd
import requests


def save_ls(ls, fp):
    with open(fp, 'w') as f:
        for i in ls:
            line = str(i) + '\n'
            f.write(line)


def download_url(tsv, save_url_path):
    df = pd.read_csv(tsv, sep='\t')
    url_ls = []
    for i in df['Download URL']:
        acc = find_acc(i)
        url = get_url(acc)
        url_ls.append(url)
    url_ls = list(set(url_ls))
    save_ls(url_ls, save_url_path)
    print(f'{len(url_ls)} zip files will be downloaded.')
    yield from url_ls


def find_acc(str):
    pattern = re.compile(r'GC[F|A]_\d*.\d')
    m = pattern.findall(str)
    if len(m) == 1:
        return m[0]
    else:
        print(str)


def get_url(acc):
    url = f'https://antismash-db.secondarymetabolites.org/output/{acc}/{acc}.zip'
    return url


def download(url, path):
    print(f'{url} is downloading')
    response = requests.get(url, stream=True)
    handle = open(os.path.join(path, os.path.basename(url)), "wb")
    # for chunk in response.iter_content(chunk_size=1024):
    #     if chunk:  # filter out keep-alive new chunks
    handle.write(response.content)
    handle.close()
    print(f'{os.path.basename(url)} download successfully!')


if __name__ == '__main__':
    pool = Pool(8)
    tsv, target_path = sys.argv[1], sys.argv[2]
    download_url = download_url(tsv, os.path.join(target_path, 'url_list.txt'))
    f = open(os.path.join(target_path, 'fail_download_list.txt'), 'w')
    for url in download_url:
        try:
            pool.apply_async(download, args=(url, target_path))
        except:
            print(str(url), file=f)
            print(f'{url} failed.')
    pool.close()
    pool.join()
    f.close()
