#! /usr/bin/env python
import argparse
import sys
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def excel2fasta(excel, seq_name_col, seq_col, fasta):
    seq_records = []
    target_col = [seq_name_col, seq_col]
    df = pd.read_excel(excel, dtype=str)[target_col]
    for index, row in df.iterrows():
        name = row[seq_name_col].strip()
        sequence = row[seq_col].strip()
        seq_records.append(SeqRecord(Seq(sequence), id=name, description=''))
    SeqIO.write(seq_records, fasta, 'fasta')


def tsv2fasta(tsv, seq_name_col, seq_col, fasta):
    seq_records = []
    target_col = [seq_name_col, seq_col]
    df = pd.read_csv(tsv, sep='\t', dtype=str)[target_col]
    for index, row in df.iterrows():
        name = row[seq_name_col].strip()
        seq = row[seq_col].replace(" ", "").strip()
        seq_record = SeqRecord(Seq(seq), id=name, description='')
        seq_records.append(seq_record)
    SeqIO.write(seq_records, fasta, "fasta")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        # prog='table2fasta',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Convert sequences stored in tsv or excel files to fasta format.")
    parser.add_argument('-i', '--input_table', required=True, type=str, help="Input tsv or excel")
    parser.add_argument('-n', '--name_col', required=True, type=str, help="Column name as sequence names in table")
    parser.add_argument('-s', '--seq_col', required=True, type=str, help="Column name as sequence in table")
    parser.add_argument('-o', '--output_fasta', required=True, type=str, help="Output fasta")
    parser.add_argument('-t', '--table_type', required=True, type=str, choices=["excel", "tsv"],
                        help="Type of input table")
    args = parser.parse_args()
    table_type = args.table_type
    if table_type == "excel":
        excel2fasta(args.input_table, args.name_col, args.seq_col, args.output_fasta)
    if table_type == "tsv":
        tsv2fasta(args.input_table, args.name_col, args.seq_col, args.output_fasta)

