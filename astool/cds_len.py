from collections import namedtuple
from astool.utils import gen_dataframe, save_dataframe2tsv, check_gbk_suffix
from astool.antismash_utils import AntismashRegionGBKParser


def get_cds_records(gbk_dir, verbose):
    # process input.
    gbk_dir_ls = []
    if check_gbk_suffix(gbk_dir):
        gbk_dir_ls.append(gbk_dir)
    else:
        with open(gbk_dir, 'r') as f:
            gbk_dir_ls = f.read().splitlines()
    gbk_dir_ls = list(set(gbk_dir_ls))

    # collect seq length information.
    for gbk in gbk_dir_ls:
        if verbose:
            print(f'{gbk} is being processed.')
        try:
            gbk_file = AntismashRegionGBKParser(gbk)
            yield gbk_file.ex_cds()
        except:
            continue


def get_cds_len(cds_records):
    """Yield cds length records for each gbk file."""
    for cds_chunk in cds_records:
        cds_lenth_records = []
        for cds in cds_chunk:
            CDSLengthRecord = namedtuple('CDSLengthRecord', ['record_id', 'locus_tag', 'length'])
            # record_id, locus_tag, translation = cds
            record_id, location, locus_tag, translation = cds
            length = len(translation)
            cds_lenth_record = CDSLengthRecord(record_id, locus_tag, length)
            cds_lenth_records.append(cds_lenth_record)
        yield cds_lenth_records


def cds_len(args):
    gbk_dir, output, verbose = args.gbk_dir, args.output, args.verbose
    cds_records = get_cds_records(gbk_dir, verbose)
    # for i in cds_records:
    #     print(i)
    cds_length_records = get_cds_len(cds_records)
    cds_length_df = gen_dataframe(cds_length_records)
    save_dataframe2tsv(cds_length_df, output)
