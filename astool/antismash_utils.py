import json
import os.path
from collections import namedtuple
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class AntismashJsonParser():
    """Generate a parser for json file from antiSMASH."""

    def __init__(self, json_dir):
        """Init.

        Args:
            json_dir: Directory of json file.
        """
        self.json_dir = json_dir

    @property
    def data(self):
        """Deserialize json to a Python object.

        Returns: A Python object.
        """
        with open(self.json_dir, 'r') as f:
            data = json.load(f)
        return data

    @property
    def input_file(self):
        return self.data['input_file']

    @property
    def records(self) -> list:
        return self.data.get('records')

    @property
    def taxon(self):
        return self.data['taxon']

    @property
    def modules_generator(self):
        for record in self.records:
            modules = record.get('modules')
            if modules:
                yield modules

    @property
    def regions(self):
        region_ls = []
        for record in self.records:
            for f in record.get("features"):
                type_ = f.get("type")
                if type_ == "region":
                    region_ls.append(f)
        return region_ls

    @property
    def cand_clusters(self):
        cand_cluster_ls = []
        for record in self.records:
            for f in record.get("features"):
                type_ = f.get("type")
                if type_ == "cand_cluster":
                    cand_cluster_ls.append(f)
        return cand_cluster_ls

    def ex_smiles_from_nrps_pks(self) -> list:
        """Extract smiles of NPRs/PKs products.

        Returns: [(record_id, region_id, smiles)]

        """
        smiles_records = []
        for record in self.records:
            try:
                record_id = record['modules']['antismash.modules.nrps_pks']['record_id']
                regions = record['modules']['antismash.modules.nrps_pks']['region_predictions']
                for region_id, smiles_list in regions.items():
                    for smiles_info in smiles_list:
                        smiles = smiles_info['smiles']
                        if smiles:
                            SmilesRecord = namedtuple("SmilesRecord", ["record_id", "region_id", "smiles"])
                            smiles_record = SmilesRecord(record_id, region_id, smiles)
                            smiles_records.append(smiles_record)
            except:
                continue

        return smiles_records

    def ex_smiles_from_nrps_pks_with_input_file(self) -> list:
        """Extract smiles of NPRs/PKs products.

        Returns: [(record_id, region_id, smiles)]

        """
        smiles_records = []
        for record in self.records:
            try:
                record_id = record['modules']['antismash.modules.nrps_pks']['record_id']
                regions = record['modules']['antismash.modules.nrps_pks']['region_predictions']
                for region_id, smiles_list in regions.items():
                    for smiles_info in smiles_list:
                        smiles = smiles_info['smiles']
                        if smiles:
                            SmilesRecord = namedtuple("SmilesRecord", ["file_path", "record_id", "region_id", "smiles"])
                            smiles_record = SmilesRecord(os.path.dirname(self.json_dir), record_id, region_id, smiles)
                            smiles_records.append(smiles_record)
            except:
                continue

        return smiles_records

    def _query_cc_type(self, record_id, region_number, cc_idx):
        records = self.records
        for record in records:
            record_id_ = record.get('id')
            if record_id_ == record_id:
                areas = record.get('areas')
                if areas:  # version = 6.0.1
                    area = areas[int(region_number) - 1]
                    protoclusters = area.get('protoclusters')
                    candidates = area.get('candidates')
                    candidate = candidates[cc_idx]
                    candidate_protoclusters = candidate.get('protoclusters')
                    cc_type_ls = []
                    for i in candidate_protoclusters:
                        cc_type_ls.append(protoclusters.get(str(i)).get('product'))
                    cc_type = "+".join(cc_type_ls)
                    # else:  # version = 5.2.0

                    return cc_type

    # def get_modules_nrps_pks_region_polymer(self) -> list:
    #     """Extract NRPS polymer prediction.
    #
    #     There are several regions in each JSON file.
    #     For each region, there is a variable number of polymer prediction results.
    #     So a dictionary was created to store the data.
    #     The keys of the dictionary are the header of table below.
    #     Finally, these dictionaries are appended to a list as return.
    #
    #     | region_id   | region_number | sc_number | polymer                      | cc_type        |
    #     | ----------- | ------------- | --------- | -----------------------------|----------------|
    #     | NC_012962.1 | 1             | 1         | (ile) + (ohmal) + (val - val)| NRPS           |
    #     | NC_012962.1 | 1             | 2         | (ile) + (ohmal) + (val - val)| NRPS+T1PKS     |
    #
    #     Returns: A list of polymer information.
    #     """
    #     polymer_list = []
    #     for modules in self.modules_generator:
    #         antismash_modules_nrps_pks = modules.get('antismash.modules.nrps_pks')
    #         if antismash_modules_nrps_pks:
    #             record_id = antismash_modules_nrps_pks.get('record_id')
    #             region_predictions = antismash_modules_nrps_pks.get('region_predictions')
    #             if record_id and region_predictions:
    #                 for region_number, cc_set in region_predictions.items():
    #                     for idx, cc in enumerate(cc_set):
    #                         # print(record_id, region_number, idx)
    #                         sc_number = cc.get('sc_number')
    #                         polymer = cc.get('polymer')
    #                         cc_type = self._query_cc_type(record_id=record_id,
    #                                                       region_number=region_number,
    #                                                       cc_idx=idx)
    #                         if polymer:
    #                             polymer_info = {
    #                                 "record_id": record_id,
    #                                 "region_number": region_number,
    #                                 "sc_number": sc_number,
    #                                 "polymer": polymer,
    #                                 "cc_type": cc_type
    #                             }
    #                             polymer_list.append(polymer_info)
    #     return polymer_list


class AntismashRegionGBKParser():
    """Generate a parser for gbk file from antiSMASH."""

    def __init__(self, gbk_dir):
        """Init.

        Args:
            gbk_dir: directory of json file.
        """
        self.gbk_dir = gbk_dir

    def _get_target_feature(self, target_type):
        for seq_record in SeqIO.parse(self.gbk_dir, 'genbank'):
            for feature in seq_record.features:
                if feature.type == target_type:
                    yield feature

    @property
    def bgc_length(self):
        """Return the length of BGC."""
        for seq_record in SeqIO.parse(self.gbk_dir, 'genbank'):
            return len(seq_record.seq)

    @property
    def bgc_type(self):
        """Return the type of BGC region."""
        bgc_type_ls = []
        for seq_record in SeqIO.parse(self.gbk_dir, 'genbank'):
            for feature in seq_record.features:
                if feature.type == 'region':
                    bgc_type_ls += feature.qualifiers['product']
        bgc_type = '+'.join(bgc_type_ls)
        return bgc_type

    @property
    def region_completeness(self):
        """Return the completeness of BGC."""
        for seq_record in SeqIO.parse(self.gbk_dir, 'genbank'):
            for feature in seq_record.features:
                if feature.type == 'region':
                    contig_edge = feature.qualifiers["contig_edge"][0]
                    if len(contig_edge) > 1:
                        return contig_edge

    @property
    def contig_seq(self):
        """Return the contig where the BGC region located."""
        for seq_record in SeqIO.parse(self.gbk_dir, 'genbank'):
            return seq_record.seq

    @property
    def contig_locus(self):
        """Return the contig where the BGC region located."""
        for seq_record in SeqIO.parse(self.gbk_dir, 'genbank'):
            return seq_record.id

    def get_nrps_pks_domains(self):
        """Return aSDomain of nrps and pks domains."""
        for seq_record in SeqIO.parse(self.gbk_dir, 'genbank'):
            for feature in seq_record.features:
                if feature.type == 'aSDomain' and feature.qualifiers.get('aSTool') == ["nrps_pks_domains"]:
                    yield feature

    def cand_cluster(self):
        for seq_record in SeqIO.parse(self.gbk_dir, 'genbank'):
            for feature in seq_record.features:
                if feature.type == 'cand_cluster':
                    CandCluster = namedtuple('CandCluster', ['location', 'smiles', 'candidate_cluster_number',
                                                             'contig_edge', 'kind', 'product'])
                    location = feature.location
                    try:
                        smiles: list = feature.qualifiers['SMILES']
                    except:
                        smiles = []
                    if smiles:
                        smiles = [''.join(smi.split()) for smi in smiles]
                    candidate_cluster_number: list = feature.qualifiers['candidate_cluster_number']
                    contig_edge: list = feature.qualifiers['contig_edge']
                    kind: list = feature.qualifiers['kind']
                    product: list = feature.qualifiers['product']
                    cand_cluster = CandCluster(location, smiles, candidate_cluster_number, contig_edge, kind, product)
                    yield cand_cluster

    def ex_smiles(self):
        """Extract all smiles from a region gbk file."""
        smiles_ls = []
        for seq_record in SeqIO.parse(self.gbk_dir, 'genbank'):
            for feature in seq_record.features:
                if feature.type == 'cand_cluster':
                    try:
                        smiles_ls += feature.qualifiers['SMILES']
                    except:
                        continue
        return smiles_ls

    def ex_domain(self):
        domain_ls = []
        for seq_record in SeqIO.parse(self.gbk_dir, 'genbank'):
            for feature in seq_record.features:
                if feature.type == 'aSDomain':
                    domain_ls += feature.qualifiers['aSDomain']
        return domain_ls

    def ex_cds(self) -> list:
        """Extract CDS sequence from the region gbk file.

        Returns: [(record_id, locus_tag, translation)]

        """
        cds_records = []
        for seq_record in SeqIO.parse(self.gbk_dir, 'genbank'):
            record_id = seq_record.id
            for feature in seq_record.features:
                if feature.type == 'CDS':
                    CDSRecord = namedtuple("CDSRecord", ["record_id", "location", "locus_tag", "translation"])
                    location = feature.location
                    locus_tag = feature.qualifiers['locus_tag'][0]
                    translation = feature.qualifiers['translation'][0]
                    cds_record = CDSRecord(record_id, location, locus_tag, translation)
                    cds_records.append(cds_record)
        return cds_records

    def ex_cds_all(self):
        cds_records = []
        for seq_record in SeqIO.parse(self.gbk_dir, 'genbank'):
            record_id = seq_record.id
            for feature in seq_record.features:
                if feature.type == 'CDS':
                    CDSRecord = namedtuple("CDSRecord", ["record_id", "location", "locus_tag", "translation",
                                                         "gene_functions", "gene_kind", "sec_met_domain"])
                    location = feature.location
                    qualifiers = feature.qualifiers
                    locus_tag = qualifiers.get('locus_tag')[0]
                    translation = qualifiers.get('translation')[0]
                    gene_functions = qualifiers.get('gene_functions')
                    gene_kind = qualifiers.get('gene_kind')
                    sec_met_domain = qualifiers.get('sec_met_domain')
                    cds_record = CDSRecord(record_id, location, locus_tag, translation,
                                           gene_functions, gene_kind, sec_met_domain)
                    cds_records.append(cds_record)
        return cds_records

    def get_cds_info_ls(self):
        seq_feat_record_ls = []
        for seq_record in SeqIO.parse(self.gbk_dir, 'genbank'):
            for feature in seq_record.features:
                if feature.type == 'CDS':
                    CDSInfo = namedtuple("CDSInfo", ['location', 'seq_record'])
                    location = feature.location
                    translation = feature.qualifiers['translation']
                    qualifiers_keys = feature.qualifiers.keys()
                    seq_id = None
                    if 'protein_id' in qualifiers_keys:
                        seq_id = feature.qualifiers['protein_id']
                    elif 'locus_tag' in qualifiers_keys:
                        seq_id = feature.qualifiers['locus_tag']
                    elif 'gene' in qualifiers_keys:
                        seq_id = feature.qualifiers['gene']
                    try:
                        gene_kind = feature.qualifiers['gene_kind']
                    except:
                        gene_kind = ['other']
                    seq_record = SeqRecord(Seq(translation[0]), id=seq_id[0], description=gene_kind[0])
                    cds_info = CDSInfo(location, seq_record)
                    seq_feat_record_ls.append(cds_info)
        return seq_feat_record_ls

    def get_nrps_pks_monomers(self):
        """Get monomers from each region.

        Returns: A dictionary of monomers items.

        """
        monomers_dict = {}
        for f in self._get_target_feature("aSDomain"):
            qualifiers = f.qualifiers
            specificity = qualifiers.get("specificity")
            if specificity:
                locus_tag = qualifiers.get("locus_tag")[0]
                for i in specificity:
                    if i.startswith("consensus"):
                        consensus = i.split(": ")[1]
                        if not monomers_dict.get(locus_tag):
                            monomers_dict[locus_tag] = [consensus]
                        else:
                            monomers_dict[locus_tag].append(consensus)
        return monomers_dict

    def ex_a_domain(self):
        """Extract sequences of A domains from each region.

        Returns:

        """
        domain_records = []
        for f in self._get_target_feature("aSDomain"):
            qualifiers = f.qualifiers
            if qualifiers.get("aSDomain") == ["AMP-binding"]:
                domain_seq = qualifiers.get("translation")[0]
                domain_id = qualifiers.get("domain_id")[0]
                consensus = "consensus: "
                consensus_ls = []
                for i in qualifiers.get("specificity"):
                    if i.startswith("consensus"):
                        consensus_ls.append(i.split(": ")[1])
                        consensus += '-'.join(consensus_ls)
                domain = SeqRecord(Seq(domain_seq), id=domain_id, description=consensus)
                domain_records.append(domain)
        return domain_records

    def ex_lanthipeptide(self):
        """Extract leader sequence and core sequence of lanthipeptide."""
        leader_seq, core_seq = [], []
        for seq_record in SeqIO.parse(self.gbk_dir, 'genbank'):
            for feature in seq_record.features:
                if feature.type == 'CDS_motif':
                    qualifiers = feature.qualifiers
                    if qualifiers.get("core_sequence") and qualifiers.get("leader_sequence"):
                        leader_seq += qualifiers.get("leader_sequence")
                        core_seq += qualifiers.get("core_sequence")
        return leader_seq, core_seq


if __name__ == '__main__':
    json_dir = "/Users/zhouzhenyi/Documents/github/astool/test_data/antiSMASH/GCA_003204095.1_ASM320409v1_genomic/GCA_003204095.1_ASM320409v1_genomic.json"
    json_file = AntismashJsonParser(json_dir)
    print(json_file.get_modules_nrps_pks_region_polymer())
