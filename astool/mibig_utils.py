import json
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


class MibigJsonParser():
    """Generate a parser for json file from MIBiG."""

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
    def mibig_accession(self):
        return self.data['cluster']['mibig_accession']

    @property
    def ncbi_tax_id(self):
        return self.data['cluster']['ncbi_tax_id']

    @property
    def organism_name(self):
        return self.data['cluster']['organism_name']

    @property
    def compounds_smiles(self):
        smiles_ls = []
        try:
            compounds = self.data['cluster']['compounds']
            for compound in compounds:
                smiles_ls.append(compound['chem_struct'])
        except:
            return smiles_ls
        return smiles_ls

    @property
    def completeness(self):
        return self.data['cluster']['loci']['completeness']

    @property
    def biosyn_class(self):
        return '+'.join(self.data['cluster']['biosyn_class'])


class MibigRegionGBKParser():
    """Generate a parser for gbk file from MIBiG."""

    def __init__(self, gbk_dir):
        """Init.

        Args:
            gbk_dir: directory of json file.
        """
        self.gbk_dir = gbk_dir

    def cds_records(self):
        seq_records = []
        for seq_record in SeqIO.parse(self.gbk_dir, 'genbank'):
            for feature in seq_record.features:
                if feature.type == 'CDS':
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
                    cds_record = SeqRecord(Seq(translation[0]), id=seq_id[0], description=gene_kind[0])
                    seq_records.append(cds_record)
        return seq_records
