import os
import sys
from zipfile import ZipFile
from fnmatch import fnmatch


def ex_json(zip_fp, target_path):
    # json_list = []
    with ZipFile(zip_fp) as zip_file:
        namelist = zip_file.namelist()
        for f in namelist:
            if f.endswith('json'):
                # json_list.append(f)
                zip_file.extract(f, target_path)
    # return json_list


def rename(fp, new_name):
    os.rename(fp, new_name)


def ex_gbk(zip_fp, target_path):
    # gbk_list = []
    with ZipFile(zip_fp) as zip_file:
        namelist = zip_file.namelist()
        for f in namelist:
            if fnmatch(f, '*region*.gbk'):
                new_name = os.path.basename(zip_fp).strip('.zip') + f'.{f}'
                # gbk_list.append(new_name)
                zip_file.extract(f, target_path)
                os.rename(os.path.join(target_path, f), os.path.join(target_path, new_name))

    # return gbk_list


def get_zip_fn_ls(folder_path):
    zip_file_list = os.listdir(folder_path)
    for f in zip_file_list:
        if fnmatch(f, 'GC?_*.zip'):
            continue
        else:
            zip_file_list.remove(f)
    print(f'共有zip文件: {len(zip_file_list)}')  # 25790
    yield from zip_file_list


if __name__ == '__main__':
    json_target_path, gbk_target_path, folder_path = sys.argv[1], sys.argv[2], sys.argv[3]
    zip_fn_ls = get_zip_fn_ls(folder_path)
    # gbk_list, json_list = [], []
    for fn in zip_fn_ls:
        try:
            ex_gbk(os.path.join(folder_path, fn), gbk_target_path)
            ex_json(os.path.join(folder_path, fn), json_target_path)
        except:
            print(f'{fn} 处理失败')
    # print(f'gbk列表长度{len(gbk_list)}，gbk集合长度{len(set(gbk_list))}')
    # print(f'json列表长度{len(json_list)}，json集合长度{len(set(json_list))}')
