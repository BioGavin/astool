#!/usr/bin/env python

import os
import sys
import zipfile
from multiprocessing import Pool

import requests


def read_url(txt):
    with open(txt, 'r') as f:
        url_ls = f.read().splitlines()
        yield from url_ls


def get_filename(url_iter):
    all_zip = []
    for url in url_iter:
        fn = url.split('/')[-1]
        all_zip.append(fn)
    return all_zip


def get_downloaded_zip(path):
    file_ls = os.listdir(path)
    downloaded_file_ls = []
    for fn in file_ls:
        fp = os.path.join(path, fn)
        if zipfile.is_zipfile(fp):
            downloaded_file_ls.append(fn)
    return downloaded_file_ls


def get_undownloaded_zip(all_zip, downloaded_zip):
    difference = list(set(all_zip).difference(set(downloaded_zip)))  # all_zip有而downloaded_zip没有的元素
    return difference


def get_url(fn):
    acc = fn[:-4]
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
    txt, download_path = sys.argv[1], sys.argv[2]
    url_iter = read_url(txt)
    all_zip = get_filename(url_iter)
    downloaded_zip = get_downloaded_zip(download_path)
    undownloaded_zip = get_undownloaded_zip(all_zip, downloaded_zip)
    while True:
        print(f'{len(undownloaded_zip)} zip files will be downloaded...')
        pool = Pool(8)
        for zip in undownloaded_zip:
            url = get_url(zip)
            try:
                pool.apply_async(download, args=(url, download_path))
            except:
                print(f'{url} failed.')
                continue
        pool.close()
        pool.join()
        downloaded_zip = get_downloaded_zip(download_path)
        undownloaded_zip = get_undownloaded_zip(all_zip, downloaded_zip)
        print(undownloaded_zip)
        if not undownloaded_zip:
            break
