from collections import namedtuple
from astool.antismash_utils import AntismashJsonParser
from astool.mibig_utils import MibigJsonParser
from astool.utils import check_json, gen_dataframe, save_dataframe2tsv


def get_json_dir_ls(input):
    """Process input."""
    json_dir_ls = []
    if check_json(input):
        json_dir_ls.append(input)
    else:
        with open(input, 'r') as f:
            json_dir_ls = f.read().splitlines()
    json_dir_ls = list(set(json_dir_ls))
    return json_dir_ls


def get_antismash_smiles_records(json_dir, verbose):
    json_dir_ls = get_json_dir_ls(json_dir)

    # collect smiles information.
    for json in json_dir_ls:
        if verbose:
            print(f'{json} is being processed.')
        try:
            json_file = AntismashJsonParser(json)
            yield json_file.ex_smiles_from_nrps_pks_with_input_file()
        except:
            continue


def get_mibig_smiles_records(json_dir, verbose):
    json_dir_ls = get_json_dir_ls(json_dir)

    # collect smiles information.
    for json in json_dir_ls:
        if verbose:
            print(f'{json} is being processed.')
        try:
            json_file = MibigJsonParser(json)
            if json_file.compounds_smiles:
                SmilesRecord = namedtuple("SmilesRecord",
                                          ["mibig_accession", "taxon_id", "biosyn_class", "completeness", "smiles"])
                smiles_records = []
                for smiles in json_file.compounds_smiles:
                    smiles_record = SmilesRecord(json_file.mibig_accession, json_file.ncbi_tax_id,
                                                 json_file.biosyn_class, json_file.completeness, smiles)
                    smiles_records.append(smiles_record)
                yield smiles_records
        except:
            continue


def ex_smiles(args):
    global smiles_records
    json_dir, output, type, verbose = args.json_dir, args.output, args.type, args.verbose
    if type == 'antismash':
        smiles_records = get_antismash_smiles_records(json_dir, verbose)
    elif type == 'mibig':
        smiles_records = get_mibig_smiles_records(json_dir, verbose)
    smiles_df = gen_dataframe(smiles_records)
    save_dataframe2tsv(smiles_df, output)
