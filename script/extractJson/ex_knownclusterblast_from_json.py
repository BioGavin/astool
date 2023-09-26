#! /usr/bin/env python3

# Usage: ex_knownclusterblast_from_json.py json_dir_list.txt knownclusterblast_hits_number_info.tsv knownclusterblast_hits_detailed_info.tsv

import os.path
import sys
import pandas as pd
from astool.antismash_utils import AntismashJsonParser
from astool.utils import save_dataframe2tsv, get_json_dir_ls


def get_knowncluster(json_dir):
    json_file = AntismashJsonParser(json_dir)
    if json_file.records:
        for record in json_file.records:
            modules = record.get("modules")
            antismash_modules_clusterblast = modules.get("antismash.modules.clusterblast")
            if antismash_modules_clusterblast:
                knowncluster = antismash_modules_clusterblast.get("knowncluster")
                yield knowncluster


if __name__ == '__main__':
    json_input, tsv_output1, tsv_output2 = sys.argv[1:4]

    knowncluster_metadata_ls = []
    knowncluster_detailed_information_ls = []

    for json_dir in get_json_dir_ls(json_input):
        json_fname = os.path.basename(json_dir)
        for kc in get_knowncluster(json_dir):
            record_id = kc.get("record_id")
            for res in kc.get("results"):
                region_number = res.get("region_number")
                total_hits = res.get("total_hits")
                data = {
                    "file_name": json_fname,
                    "record_id": record_id,
                    "region_number": region_number,
                    "total_hits": total_hits
                }
                knowncluster_metadata_ls.append(data)
                ranking = res.get("ranking")
                if total_hits > 0:
                    for rank in ranking:
                        accession, description, cluster_type = rank[0].get("accession"), \
                                                               rank[0].get("description"), \
                                                               rank[0].get("cluster_type")
                        known_bgc_gen_len = len(rank[0].get("proteins"))
                        pairings = rank[1].get("pairings")
                        blast_score = rank[1].get("blast_score")
                        hit_known_gene = []
                        for pair in pairings:
                            hit_known_gene_name = pair[2].get("name")
                            hit_known_gene.append(hit_known_gene_name)
                        hit_known_gene = list(set(hit_known_gene))
                        proportion = len(hit_known_gene) / known_bgc_gen_len * 100
                        detailed_data = {
                            "file_name": json_fname,
                            "record_id": record_id,
                            "region_number": region_number,
                            "MIBiG_accession": accession,
                            "MIBiG_description": description,
                            "MIBiG_cluster_type": cluster_type,
                            "blast_score": blast_score,
                            "MIBiG_similarity_precise_proportion": proportion,
                            "MIBiG_similarity_rough_proportion": int(proportion)
                        }
                        knowncluster_detailed_information_ls.append(detailed_data)

    save_dataframe2tsv(pd.DataFrame(knowncluster_metadata_ls), tsv_output1)
    save_dataframe2tsv(pd.DataFrame(knowncluster_detailed_information_ls), tsv_output2)
