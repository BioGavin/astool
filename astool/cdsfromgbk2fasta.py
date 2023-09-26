from Bio import SeqIO
from astool.mibig_utils import MibigRegionGBKParser

def cdsfromgbk2fasta(args):
    gbk_dir, fasta_dir = args.gbk_dir, args.fasta_dir
    gbk_file = MibigRegionGBKParser(gbk_dir)
    seq_records = gbk_file.cds_records()
    SeqIO.write(seq_records, fasta_dir, "fasta")