#!/usr/bin/env python3
#Usage: python3 ex_polymer_from_json.py json.list output.tsv

# """Extract NRPS polymer prediction.
# Applicable antiSMASH version 5.2.0
# There are several regions in each JSON file.
# For each region, there is a variable number of polymer prediction results.
# So a dictionary was created to store the data.
# The keys of the dictionary are the header of table below.
# Finally, these dictionaries are appended to a list as return.
#
# | json_file    | record_id | region_number | region_type | cc_number | cc_type        | polymer                      |
# | ------------ | --------- | ------------- | ----------- | --------- | -------------- | ---------------------------- |
# | GCA0000.json | NC_0129.1 | 1             | NRPS+T1PKS  | 1         | NRPS           | (ile) + (ohmal) + (val - val)|
# | GCA0000.json | NC_0129.1 | 1             | NRPS+T1PKS  | 2         | NRPS+T1PKS     | (ile) + (ohmal) + (val - val)|
#
# """
import os.path
import sys
import pandas as pd
from astool.antismash_utils import AntismashJsonParser
from astool.utils import get_json_dir_ls, save_dataframe2tsv, check_json


def get_candidate_cluster_number_and_cc_type_dict(record):
    data_dict = {}
    for f in record.get("features"):
        type_ = f.get("type")
        if type_ == "cand_cluster":
            qualifiers = f.get("qualifiers")
            candidate_cluster_number = qualifiers.get("candidate_cluster_number")[0]
            cc_type = "+".join(qualifiers.get("product"))
            data_dict[candidate_cluster_number] = cc_type
    return data_dict


def get_polymer(record, region_number, cc_number: str):
    modules = record.get("modules")
    nrps_pks_modules = modules.get("antismash.modules.nrps_pks")
    if nrps_pks_modules:
        region_predictions = nrps_pks_modules.get("region_predictions")
        cc_data = region_predictions.get(region_number)
        if cc_data:
            for cc in cc_data:  # 这里可能有问题
                if cc.get("sc_number") == int(cc_number):
                    polymer = cc.get("polymer")
                    return polymer


def main(json):
    data_ls = []
    json_file = json.json_dir
    records = json.records
    if records:
        for record in records:
            cc_number_and_type_dict = get_candidate_cluster_number_and_cc_type_dict(record)
            record_id = record.get("id")
            features = record.get("features")
            for feature in features:
                if feature.get("type") == "region":
                    qualifiers = feature.get("qualifiers")
                    region_number = qualifiers.get("region_number")[0]
                    region_type = "+".join(qualifiers.get("product"))
                    candidate_cluster_numbers = qualifiers.get("candidate_cluster_numbers")
                    for cc_number in candidate_cluster_numbers:
                        cc_type = cc_number_and_type_dict.get(cc_number)
                        if region_number and cc_number:
                            polymer = get_polymer(record, region_number, cc_number)
                            if polymer:
                                data = {
                                    "json_file": json_file,
                                    "record_id": record_id,
                                    "region_number": region_number,
                                    "region_type": region_type,
                                    "cc_number": cc_number,
                                    "cc_type": cc_type,
                                    "polymer": polymer
                                }

                                data_ls.append(data)
        return data_ls


if __name__ == '__main__':
    json_input, tsv_output = sys.argv[1], sys.argv[2]

    data_ls = []
    failed_jobs = []

    for json_dir in get_json_dir_ls(json_input):
        if check_json(json_dir):
            json_file = AntismashJsonParser(json_dir)
            data_ls += main(json_file)
        else:
            failed_jobs.append(json_dir + "\n")

    with open("log.txt", "w") as f:
        f.writelines(failed_jobs)
    save_dataframe2tsv(pd.DataFrame(data_ls), tsv_output)
