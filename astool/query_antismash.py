# Query antiSMASH file (GBK or JSON) contents by keyword or ID
import json
from astool.antismash_utils import AntismashJsonParser


class AntismashJsonQueryer(AntismashJsonParser):

    def __init__(self, json_dir):
        super().__init__(json_dir)


    def query_cc_type(self, record_id, region_number, sc_number_idx):
        records = self.records
        for record in records:
            record_id_ = record.get('id')
            if record_id_ == record_id:
                areas = record.get('areas')
                area = areas[int(region_number)-1]
                protoclusters = area.get('protoclusters')
                candidates = area.get('candidates')
                candidate = candidates[sc_number_idx]
                candidate_protoclusters = candidate.get('protoclusters')
                cc_type_ls = []
                for i in candidate_protoclusters:
                    cc_type_ls.append(protoclusters.get(str(i)).get('product'))
                cc_type = "+".join(cc_type_ls)
                return cc_type









if __name__ == '__main__':
    # test()
    json_dir = '/Users/zhouzhenyi/Documents/SCIProject/BiLinker/astool/test_data/antiSMASH/GCF_000196475.1_ASM19647v1_genomic/GCF_000196475.1_ASM19647v1_genomic.json'
    queryer = AntismashJsonQueryer(json_dir)
    region_id = 'NC_012962.1'
    print(queryer.query_cc_type(region_id, '1', '1'))
