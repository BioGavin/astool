import argparse
import os

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from astool.antismash_utils import AntismashRegionGBKParser
from Bio.Seq import Seq

from astool.utils import check_gbk_suffix, save_dataframe2tsv


def get_gbk_dir(input_file):
    assert os.path.isfile(input_file), "The input is not a file."
    gbk_list = []
    if check_gbk_suffix(input_file):
        gbk_list.append(input_file)
    else:
        with open(input_file, 'r') as f:
            antismash_folder_list = f.read().splitlines()
            for antismash_folder in antismash_folder_list:
                region_gbk = [i for i in os.listdir(antismash_folder) if "region" in i and i.endswith("gbk")]
                region_gbk_dir = [os.path.join(antismash_folder, i) for i in region_gbk]
                gbk_list += region_gbk_dir
    return gbk_list


def ex_contig_and_cds_for_single(gbk_file):
    gbk = AntismashRegionGBKParser(gbk_file)
    contig_seq = gbk.contig_seq
    contig_locus = gbk.contig_locus
    cds_records = gbk.ex_cds_all()
    return contig_locus, contig_seq, cds_records


def write_cds_fna(cds_records, contig_seq, output_fna):
    records = []
    for cds in cds_records:
        if cds.gene_functions:
            description = " | ".join(cds.gene_functions)
        else:
            description = ""
        location = cds.location
        record = SeqRecord(seq=Seq(contig_seq[location.start: location.end]), id=cds.locus_tag, description=description)
        records.append(record)
    SeqIO.write(records, output_fna, "fasta")


def write_cds_faa(cds_records, output_faa):
    records = []
    for cds in cds_records:
        if cds.gene_functions:
            description = " | ".join(cds.gene_functions)
        else:
            description = ""
        record = SeqRecord(seq=Seq(cds.translation), id=cds.locus_tag, description=description)
        records.append(record)
    SeqIO.write(records, output_faa, "fasta")


def gen_cds_df(assembly_accession, contig_locus, cds_records):
    assembly_accession_ls = []
    contig_locus_ls = []
    cds_location_ls = []
    cds_locus_tag_ls = []
    gene_kind_ls = []
    gene_functions_ls = []
    sec_met_domain_ls = []
    for cds in cds_records:
        assembly_accession_ls.append(assembly_accession)
        contig_locus_ls.append(contig_locus)
        cds_location_ls.append(cds.location)
        cds_locus_tag_ls.append(cds.locus_tag)
        if cds.gene_kind:
            gene_kind_ls.append(" | ".join(cds.gene_kind))
        else:
            gene_kind_ls.append("")
        if cds.gene_functions:
            gene_functions_ls.append(" | ".join(cds.gene_functions))
        else:
            gene_functions_ls.append("")
        if cds.sec_met_domain:
            sec_met_domain_ls.append(" | ".join(cds.sec_met_domain))
        else:
            sec_met_domain_ls.append("")

    df = pd.DataFrame({
        'assembly_accession': assembly_accession_ls,
        'contig_locus': contig_locus_ls,
        'cds_location': cds_location_ls,
        'cds_locus_tag': cds_locus_tag_ls,
        'gene_kind': gene_kind_ls,
        'gene_functions': gene_functions_ls,
        'sec_met_domain': sec_met_domain_ls
    })
    return df


def help():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
examples:
1. python ex_CDS_from_region_gbk.py -i GG657748.1.region001.gbk -o CDS
2. python ex_CDS_from_region_gbk.py -i antismash_result_folder.list -o CDS
        '''
    )

    parser.add_argument('-i', '--input', required=True, type=str,
                        help='Input a single region GBK file or a list of antiSMASH result folder paths.')
    parser.add_argument('-o', '--output', required=True, type=str,
                        help='Output folder paths.')
    parser.add_argument('-f', '--format', required=False, action='store_true',
                        help='If the antiSMASH result folder name begins with assembly accession'
                             '(like "GCA_000000001.1" or "GCF_000000001.1"),'
                             'enable the flag to format the output file name.')

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = help()
    gbk_dir_list = get_gbk_dir(args.input)
    if len(gbk_dir_list) == 1:  # 处理单个Region的GBK文件
        contig_locus, contig_seq, cds_records = ex_contig_and_cds_for_single(gbk_dir_list[0])

        write_cds_faa(cds_records,
                      os.path.join(args.output, os.path.basename(gbk_dir_list[0]).replace("gbk", "cds.faa")))

        write_cds_fna(cds_records, contig_seq,
                      os.path.join(args.output, os.path.basename(gbk_dir_list[0]).replace("gbk", "cds.fna")))
        genome_file_name = os.path.basename(os.path.dirname(gbk_dir_list[0]))
        if args.format:
            genome_file_name = genome_file_name[:15]
        cds_df = gen_cds_df(genome_file_name, contig_locus, cds_records)
        save_dataframe2tsv(cds_df,
                           os.path.join(args.output, os.path.basename(gbk_dir_list[0]).replace("gbk", "cds.tsv")))
    else:  # 批量处理antiSMASH结果中Region的GBK文件
        faa_output_folder = os.path.join(args.output, "faa")
        fna_output_folder = os.path.join(args.output, "fna")
        tsv_output_folder = os.path.join(args.output, "tsv")

        os.makedirs(faa_output_folder, exist_ok=True)
        os.makedirs(fna_output_folder, exist_ok=True)
        os.makedirs(tsv_output_folder, exist_ok=True)

        for gbk_dir in gbk_dir_list:  # 循环处理每个Region的GBK文件
            contig_locus, contig_seq, cds_records = ex_contig_and_cds_for_single(gbk_dir)
            gbk_file_name = os.path.basename(gbk_dir)
            genome_file_name = os.path.basename(os.path.dirname(gbk_dir))
            if args.format:
                genome_file_name = genome_file_name[:15]
            write_cds_faa(cds_records,
                          os.path.join(faa_output_folder,
                                       genome_file_name + "." + gbk_file_name.replace("gbk", "cds.faa")))
            write_cds_fna(cds_records, contig_seq,
                          os.path.join(fna_output_folder,
                                       genome_file_name + "." + gbk_file_name.replace("gbk", "cds.fna")))
            cds_df = gen_cds_df(genome_file_name, contig_locus, cds_records)
            save_dataframe2tsv(cds_df,
                               os.path.join(tsv_output_folder,
                                            genome_file_name + "." + gbk_file_name.replace("gbk", "cds.tsv")))
