# Extract candidate cluster from gbk file and generate fasta file.
# Extract SMILES from each candidate cluster.


import os.path
import sys

import pandas as pd
from Bio import SeqIO

from astool.utils import check_gbk
from astool.antismash_utils import AntismashRegionGBKParser
import astool


def get_gbk_dir_ls(input):
    """Process input."""
    gbk_dir_ls = []
    if check_gbk(input):
        gbk_dir_ls.append(input)
    else:
        with open(input, 'r') as f:
            gbk_dir_ls = f.read().splitlines()
    gbk_dir_ls = list(set(gbk_dir_ls))
    yield from gbk_dir_ls


def gen_cc_fasta(cds_info_ls, cc_info, fasta_dir):
    """Generate fasta file for candidate cluster.

    Args:
        gbk_file:
        cc_info:

    Returns: A file in fasta format.

    """
    cc_location_start, cc_location_end = cc_info.location.start, cc_info.location.end
    cc_cds_records = []
    for cds_info in cds_info_ls:
        if cc_location_start <= cds_info.location.start and cds_info.location.end <= cc_location_end:
            cc_cds_records.append(cds_info.seq_record)
    SeqIO.write(cc_cds_records, fasta_dir, "fasta")


if __name__ == '__main__':
    gbk_list_fp, fasta_folder = sys.argv[1], sys.argv[2]
    gbk_dir_ls = get_gbk_dir_ls(gbk_list_fp)
    genome_id_ls = []
    gbk_id_ls = []
    smi_ls = []
    bgc_type_ls = []
    completeness_ls = []
    for gbk_dir in gbk_dir_ls:
        gbk_file = AntismashRegionGBKParser(gbk_dir)
        cds_info_ls = gbk_file.get_cds_info_ls()
        for cc_info in gbk_file.cand_cluster():
            if cc_info.smiles and cc_info.smiles[0]:
                # 为每个candidate_cluster生成fasta文件
                fasta_fn = os.path.basename(gbk_dir).strip('.gbk') + '.cc' + \
                           cc_info.candidate_cluster_number[0] + '.fasta'
                fasta_dir = os.path.join(fasta_folder, fasta_fn)
                gen_cc_fasta(cds_info_ls, cc_info, fasta_dir)
                # 收集metadata信息
                genome_id_ls.append(os.path.basename(gbk_dir)[:15])
                gbk_id_ls.append(fasta_fn.strip('.fasta'))
                smi_ls.append(cc_info.smiles[0])
                bgc_type_ls.append('+'.join(cc_info.product))
                completeness_ls.append(cc_info.contig_edge[0])
                print(os.path.basename(gbk_dir)[:15],
                      fasta_fn.strip('.fasta'),
                      cc_info.smiles[0],
                      '+'.join(cc_info.product),
                      cc_info.contig_edge[0])
    pd.DataFrame({
        'genome_id': genome_id_ls,
        'gbk_id': gbk_id_ls,
        'smiles': smi_ls,
        'completeness': completeness_ls,
        'bgc_type': bgc_type_ls
    }).to_csv('smiles.tsv', index=False, sep='\t')
