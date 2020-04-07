import numpy as np

from idbd_bio_utils import NcbiTaxonomy

class SampleCompParser(object):
    def __init__(self, sample_composition_path, ncbi_tax=None,
            ctrl_taxa=None):
        if ncbi_tax is None:
            self.ncbi_tax = NcbiTaxonomy()
        else:
            self.ncbi_tax = ncbi_tax
        self._validate_ctrl_taxa(ctrl_taxa)
        self.composition_dict = self._load_sample_composition(
            sample_composition_path)
        self.total_reads = np.sum(list(self.composition_dict.values()))
        self.organism_counts = self._get_organism_counts()

    def _validate_ctrl_taxa(self, ctrl_taxa):
        if isinstance(ctrl_taxa, list):
            if not all([isinstance(taxid, int) for taxid in ctrl_taxa]):
                raise ValueError('ctrl_taxa must be int or list of ints.')
            self.ctrl_taxa = ctrl_taxa
        elif isinstance(ctrl_taxa, int):
            self.ctrl_taxa = [ctrl_taxa]
        elif ctrl_taxa is None:
            self.ctrl_taxa = [ctrl_taxa]
        else:
            raise ValueError('ctrl_taxa must be int or list of ints.')

    def _validate_taxid(self, taxid):
        if not isinstance(taxid, int):
            try:
                return int(taxid)
            except:
                raise ValueError(f'taxid must be int, received {type(taxid)}')
        else:
            return taxid

    def _load_sample_composition(self, sample_composition_path):
        with open(sample_composition_path) as infile:
            composition_dict = {}
            for line in infile:
                data = line.strip().split('\t')
                composition_dict.update({int(data[0]): int(data[1])})
        return composition_dict

    def _get_organism_counts(self):
        human_taxid = 9606
        bacteria_taxid = 2
        virus_taxid = 10239
        eukaryota_taxid = 2759
        parasite_taxids = {
            7563,
            188941,
            6029,
            5653,
            6935,
            6178,
            5794,
            6308,
            31277,
            119088,
            6199,
            85819,
            33083,
            33084,
            75966,
            41165,
            7509,
            6236,
            198624,
            33634,
            5988,
            6249,
            5738,
            1489900,
            740972,
            1485168,
            37104,
            10232
        }
        fungus_taxid = 4751

        human_count = 0
        bacteria_count = 0
        virus_count = 0
        parasite_count = 0
        fungus_count = 0
        unclassified_count = 0
        for taxid, count in self.composition_dict.items():
            if taxid in self.ctrl_taxa:
                continue
            taxid_path = self.ncbi_tax.get_path(taxid)
            if human_taxid in taxid_path:
                human_count += count
            elif fungus_taxid in taxid_path:
                fungus_count += count
            elif parasite_taxids.intersection(set(taxid_path)) != set():
                parasite_count += count
            elif bacteria_taxid in taxid_path:
                bacteria_count += count
            elif virus_taxid in taxid_path:
                virus_count += count
            else:
                unclassified_count += count
        count_dict = {
            'human': human_count,
            'bacteria': bacteria_count,
            'virus': virus_count,
            'parasite': parasite_count,
            'fungus': fungus_count,
            'unclassified': unclassified_count
        }
        return count_dict

    def get_total_reads(self):
        return self.total_reads

    def get_taxid_reads(self, taxid):
        taxid = self._validate_taxid(taxid)
        if taxid not in self.composition_dict:
            reads = 0
        else:
            reads = self.composition_dict[taxid]
        return reads

    def get_taxid_nr(self, taxid, normalizer=1e7):
        taxid = self._validate_taxid(taxid)
        if taxid not in self.composition_dict:
            nr = 0
        else:
            nr = normalizer * self.composition_dict[taxid] / self.total_reads
        return nr

    def get_genus_nr(self, genus_taxid, normalizer=1e7):
        genus_taxid = self._validate_taxid(genus_taxid)
        children = self.ncbi_tax.get_children(genus_taxid) + [genus_taxid]
        intersection = set(children).intersection(
            set(list(self.composition_dict.values())))
        reads = 0
        for taxid in intersection:
            if taxid in self.composition_dict:
                reads += self.composition_dict[taxid]
        # nr = normalizer * reads / self.total_reads
        nr = reads
        return nr

    def get_org_comp_abs(self):
        return self.organism_counts

    def get_org_comp_rel(self):
        count_dict = {}
        for org, count in self.organism_counts.items():
            count_dict.update({org: 100 * count / self.total_reads})
        return count_dict

    def get_org_comp_nr(self, normalizer=1e7):
        count_dict = {}
        for org, count in self.organism_counts.items():
            count_dict.update({org: normalizer * count / self.total_reads})
        return count_dict
