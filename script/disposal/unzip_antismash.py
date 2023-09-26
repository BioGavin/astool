import os
import zipfile


def get_unzip_ls(zip_db_pf, unzipped_db_pf):
    # 获取待解压的文件列表
    zip_db = os.listdir(zip_db_pf)
    zip_db = [f.strip('.zip') for f in zip_db]
    unzipped_db = os.listdir(unzipped_db_pf)
    difference = list(set(zip_db).difference(set(unzipped_db)))  # zip_db有而unzipped_db没有的元素
    return difference


def unzip_file(zip_fp, target_path):
    zip_file = zipfile.ZipFile(zip_fp, 'r')
    unzip_fn = os.path.basename(zip_fp).strip('.zip')
    unzip_fp = os.path.join(target_path, unzip_fn)
    zip_file.extractall(unzip_fp)


def unzip_antismash(zip_db_pf, unzipped_db_pf):
    unzip_ls = get_unzip_ls(zip_db_pf, unzipped_db_pf)
    print(f'{len(unzip_ls)}个文件待解压...')
    for f in unzip_ls:
        fn = f + '.zip'
        zip_fp = os.path.join(zip_db_pf, fn)
        try:
            unzip_file(zip_fp, unzipped_db_pf)
        except:
            print(f'{f}解压失败')


if __name__ == '__main__':
    zip_db_pf = '/Volumes/W7/antismash_db_results'
    unzipped_db_pf = '/Volumes/W7/antismash_db'
    unzip_antismash(zip_db_pf, unzipped_db_pf)
