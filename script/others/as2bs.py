#! /usr/bin/env python3

import argparse

import pandas as pd


def read_input(input_file):
    if input_file.endswith("xlsx"):
        input_df = pd.read_excel(input_file)
    elif input_file.endswith("tsv"):
        input_df = pd.read_csv(input_file, sep='\t')
    else:
        raise ValueError("The input file is not in the specified format (tsv or xlsx).")
    columns_ls = input_df.columns.to_list()
    region_num_idx = columns_ls.index("Region_Num")
    static_df = input_df[columns_ls[:region_num_idx + 1]]
    dynamic_df = input_df[columns_ls[region_num_idx + 1:]]

    return static_df, dynamic_df


def convert_bgc(product):
    """Sort BGC by its type. Uses antiSMASH annotations
    (see
    https://docs.antismash.secondarymetabolites.org/glossary/#cluster-types)
    """

    # TODO: according with current (2021-05) antiSMASH rules:
    # prodigiosin and PpyS-KS -> PKS
    # CDPS, mycosporine-like -> NRPS
    # but they should probably be kept in Others
    pks1_products = {'t1pks', 'T1PKS'}
    pksother_products = {'transatpks', 't2pks', 't3pks', 'otherks', 'hglks',
                         'transAT-PKS', 'transAT-PKS-like', 'T2PKS', 'T3PKS',
                         'PKS-like', 'hglE-KS', 'prodigiosin'}
    nrps_products = {'nrps', 'NRPS', 'NRPS-like', 'thioamide-NRP',
                     'NAPAA'}
    # Zhen-Yi Zhou(gavinchou99@126.com) added lanthidin to ripps_products
    ripps_products = {'lantipeptide', 'thiopeptide', 'bacteriocin', 'linaridin',
                      'cyanobactin', 'glycocin', 'LAP', 'lassopeptide',
                      'sactipeptide', 'bottromycin', 'head_to_tail', 'microcin',
                      'microviridin', 'proteusin', 'lanthipeptide', 'lipolanthine',
                      'RaS-RiPP', 'fungal-RiPP', 'TfuA-related', 'guanidinotides',
                      'RiPP-like', 'lanthipeptide-class-i', 'lanthipeptide-class-ii',
                      'lanthipeptide-class-iii', 'lanthipeptide-class-iv',
                      'lanthipeptide-class-v', 'ranthipeptide', 'redox-cofactor',
                      'thioamitides', 'epipeptide', 'cyclic-lactone-autoinducer',
                      'spliceotide', 'RRE-containing', 'crocagin', 'lanthidin'}
    saccharide_products = {'amglyccycl', 'oligosaccharide', 'cf_saccharide',
                           'saccharide'}
    others_products = {'acyl_amino_acids', 'arylpolyene', 'aminocoumarin',
                       'ectoine', 'butyrolactone', 'nucleoside', 'melanin',
                       'phosphoglycolipid', 'phenazine', 'phosphonate', 'other',
                       'cf_putative', 'resorcinol', 'indole', 'ladderane',
                       'PUFA', 'furan', 'hserlactone', 'fused', 'cf_fatty_acid',
                       'siderophore', 'blactam', 'fatty_acid', 'PpyS-KS', 'CDPS',
                       'betalactone', 'PBDE', 'tropodithietic-acid', 'NAGGN',
                       'halogenated', 'pyrrolidine', 'mycosporine-like'}

    # PKS_Type I
    if product in pks1_products:
        return ("PKSI")
    # PKS Other Types
    elif product in pksother_products:
        return ("PKSother")
    # NRPs
    elif product in nrps_products:
        return ("NRPS")
    # RiPPs
    elif product in ripps_products:
        return ("RiPPs")
    # Saccharides
    elif product in saccharide_products:
        return ("Saccharides")
    # Terpenes
    elif product == 'terpene':
        return ("Terpene")
    # PKS/NRP hybrids
    elif len(product.split("+")) > 1:
        # print("  Possible hybrid: (" + cluster + "): " + product)
        # cf_fatty_acid category contains a trailing empty space
        subtypes = set(s.strip() for s in product.split("+"))
        if len(subtypes - (pks1_products | pksother_products | nrps_products)) == 0:
            if len(subtypes - nrps_products) == 0:
                return ("NRPS")
            elif len(subtypes - (pks1_products | pksother_products)) == 0:
                return ("PKSother")  # pks hybrids
            else:
                return ("PKS-NRP_Hybrids")
        elif len(subtypes - ripps_products) == 0:
            return ("RiPPs")
        elif len(subtypes - saccharide_products) == 0:
            return ("Saccharides")
        else:
            return ("Others")  # other hybrid
    # Others
    elif product in others_products:
        return ("Others")
    # ??
    elif product == "":
        # No product annotation. Perhaps not analyzed by antiSMASH
        return ("Others")
    else:
        print("  Warning: unknown product '{}'".format(product))
        return ("Others")


def convert_df(df):
    catgory_dict = {"PKSI": [], "PKSother": [], "NRPS": [], "RiPPs": [],
                    "Saccharides": [], "Terpene": [], "PKS-NRP_Hybrids": [], "Others": []}
    columns_ls = df.columns
    for column in columns_ls:
        catgory_dict[convert_bgc(column.strip())].append(column)

    new_df = pd.DataFrame()
    for k, v in catgory_dict.items():
        new_df[k] = dynamic_df[v].sum(axis=1)

    return new_df


def save_df(output_file, tosave_df):
    if output_file.endswith("xlsx"):
        tosave_df.to_excel(output_file, index=False)
    else:
        tosave_df.to_csv(output_file, sep='\t', index=False)


def parse_args():
    parser = argparse.ArgumentParser(
        prog='as2bs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Convert antiSMASH-BGC-Types to BiGSPACE-BGC-Types")
    parser.add_argument('-i', '--input', required=True, type=str, help="Input tsv or excel")
    parser.add_argument('-o', '--output', required=True, type=str, help="Output tsv or excel")
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = parse_args()
    static_df, dynamic_df = read_input(args.input)
    new_df = convert_df(dynamic_df)
    out_df = pd.concat([static_df, new_df], axis=1)
    save_df(args.output, out_df)

