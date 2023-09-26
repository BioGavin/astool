import os
import re
import sys

import pandas as pd
from astool.bgc_info import BGC


def get_region_gbk_dirs(antismash_folder):
    dirs = os.listdir(antismash_folder)

    gbk_dirs = [os.path.join(antismash_folder, d) for d in dirs if "region" in d and d.endswith("gbk")]

    return gbk_dirs


if __name__ == '__main__':
    antismash_folder, output_tsv = sys.argv[1:3]
    # antismash_folder = "/Users/zhouzhenyi/Documents/github/wlabkit/astool/test_data/antiSMASH/GCF_002968995.1_ASM296899v1_genomic"
    pattern = r"(GC[AF]_\d{9}\.\d)"
    genome_id = re.findall(pattern, antismash_folder)[0]

    gbk_dirs = get_region_gbk_dirs(antismash_folder)
    gene_info_df_ls = []
    pfam_info_df_ls = []
    for gbk_dir in gbk_dirs:
        bgc = BGC(gbk_dir)
        contig_id = bgc.region_id
        gene_info_data_ls = []
        for gene_info in bgc.build_gene_info():
            gene_info_data = {}
            locus_tag = gene_info.locus_tag
            gene_functions = gene_info.gene_functions
            location = gene_info.location
            translation = gene_info.translation
            gene_info_data["genome_id"] = genome_id
            gene_info_data["contig_id"] = contig_id
            gene_info_data["gene_id"] = locus_tag
            gene_info_data["location"] = location
            gene_info_data["gene_function"] = gene_functions
            gene_info_data["AA_seq"] = translation
            gene_info_data_ls.append(gene_info_data)
        gene_info_df_ls.append(pd.DataFrame(gene_info_data_ls))

        pfam_info_data_ls = []
        for pfam_info in bgc.build_pfam_info():
            pfam_info_data = {}
            locus_tag = pfam_info.locus_tag
            pfam_hit = pfam_info.pfam_hit

            pfam_info_data["gene_id"] = locus_tag
            pfam_info_data["pfam_hit"] = pfam_hit

            pfam_info_data_ls.append(pfam_info_data)
        pfam_info_df_ls.append(pd.DataFrame(pfam_info_data_ls))

    gene_info_df = pd.concat(gene_info_df_ls)

    pfam_info_df = pd.concat(pfam_info_df_ls)
    if pfam_info_df.empty:
        pfam_info_df = pd.DataFrame(columns=["gene_id", "pfam_hits"])
    # print(gene_info_df)
    # print("---")
    # print(pfam_info_df)
    final_info_df = pd.merge(left=gene_info_df, right=pfam_info_df, how="left", left_on=["gene_id"],
                             right_on=["gene_id"])

    # gene_info_df.set_index("gene_id", inplace=True)
    # pfam_info_df.set_index("gene_id", inplace=True)
    # final_info_df = pd.concat([gene_info_df, pfam_info_df], axis=1)
    final_info_df.to_csv(output_tsv, sep='\t', index=False)
