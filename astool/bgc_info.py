from Bio import SeqIO
from astool.antismash_utils import AntismashRegionGBKParser


class BGC(AntismashRegionGBKParser):
    def __init__(self, gbk_dir):
        super().__init__(gbk_dir)
        self.content = SeqIO.parse(self.gbk_dir, 'genbank')

    @property
    def region_id(self):
        region_id = ''
        for seq_record in self.content:
            region_id = seq_record.id
        return region_id

    def build_gene_info(self):
        for cds_record in self._get_target_feature("CDS"):
            cds_gene = GeneInfo(cds_record)
            yield cds_gene

    def build_pfam_info(self):
        for pfam_record in self._get_target_feature("PFAM_domain"):
            pfam = PfamInfo(pfam_record)
            yield pfam


class GeneInfo(BGC):
    def __init__(self, cds_record):
        self.cds_record = cds_record
        self.qualifiers = cds_record.qualifiers

    @property
    def locus_tag(self):
        locus_tag = self.qualifiers.get("locus_tag")[0]
        return locus_tag

    @property
    def location(self):
        location = self.cds_record.location
        return location

    @property
    def translation(self):
        translation = self.qualifiers.get("translation")[0]
        return translation

    @property
    def gene_functions(self):
        gene_functions = self.qualifiers.get("gene_functions")
        if gene_functions:
            gene_functions = gene_functions[0]
        return gene_functions


class PfamInfo(BGC):
    def __init__(self, pfam_record):
        self.pfam_record = pfam_record
        self.qualifiers = pfam_record.qualifiers

    @property
    def locus_tag(self):
        locus_tag = self.qualifiers.get("locus_tag")[0]
        return locus_tag

    @property
    def db_xref(self):
        db_xref = self.qualifiers.get("db_xref")[0]
        return db_xref

    @property
    def description(self):
        description = self.qualifiers.get("description")[0]
        return description

    @property
    def evalue(self):
        evalue = self.qualifiers.get("evalue")[0]
        return evalue

    @property
    def score(self):
        score = self.qualifiers.get("score")[0]
        return score

    @property
    def protein_start(self):
        protein_start = self.qualifiers.get("protein_start")[0]
        return protein_start

    @property
    def protein_end(self):
        protein_end = self.qualifiers.get("protein_end")[0]
        return protein_end

    @property
    def pfam_hit(self):
        pfam_hit = f"{self.db_xref} ({self.description}): [{self.protein_start}:{self.protein_end}] (score: {self.score}, e-value: {self.evalue})"
        return pfam_hit


if __name__ == '__main__':
    test_gbk = "/Users/zhouzhenyi/Downloads/tmp/antismash/GCA_000525635.1_ASM52563v1_genomic/CP007155.1.region048.gbk"
    bgc = BGC(test_gbk)

    # for gene_info in bgc.build_gene_info():
    #     print(gene_info.gene_functions)
    for pfam_info in bgc.build_pfam_info():
        print(pfam_info.locus_tag, pfam_info.pfam_hit)

    # print(pfam_info.locus_tag)
    # print(pfam_info.pfam_hits)
