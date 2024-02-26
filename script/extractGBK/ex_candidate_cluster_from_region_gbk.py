import argparse
import os.path
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def read_region_gbk(gbk):
    for seq_record in SeqIO.parse(gbk, 'genbank'):
        return seq_record


def ex_cc(gbk, cc_number):
    gbk_basename = os.path.basename(gbk)

    seq_record = None
    for record in SeqIO.parse(gbk, "genbank"):
        seq_record = record
        break
    seq = seq_record.seq
    seq_id = seq_record.id
    features = seq_record.features
    for feature in features:  # 依次读取每个 Feature
        if feature.type == 'cand_cluster':
            location = feature.location
            # print(feature.qualifiers)
            candidate_cluster_number = feature.qualifiers.get('candidate_cluster_number')
            kind = " ".join(feature.qualifiers.get('kind'))
            if cc_number in candidate_cluster_number:
                target_seq = seq[location.start:location.end]
                target_seq_record = SeqRecord(Seq(target_seq), id=seq_id,
                                              description=f"candidate_cluster {cc_number} form {gbk_basename}")
                return target_seq_record


def save_target_seq(seq, fna):
    SeqIO.write(seq, fna, 'fasta')


def help():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
examples:
1. python ex_candidate_cluster_from_region_gbk.py -i Genome000001.region001.gbk -n 1 -o cc.fna
        '''
    )

    parser.add_argument('-i', '--input', required=True, type=str,
                        help='Input a single region GBK file.')
    parser.add_argument('-o', '--output', required=True, type=str,
                        help='Output candidate cluster sequence in FASTA format.')
    parser.add_argument('-n', '--cc_number', required=True, type=str)

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = help()

    target_seq_record = ex_cc(args.input, args.cc_number)
    save_target_seq(target_seq_record, args.output)
