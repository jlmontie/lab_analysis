import sys

import pandas as pd
import numpy as np
from tqdm import tqdm

from idbd_bio_utils import NcbiTaxonomy
sys.path.insert(0, '../scripts')
from sample_composition_utils import SampleCompParser


def get_quant(sample_info, panel_orgs, coverages):
    # matrix = np.zeros(len(seq_sple_ls), len(panel_orgs))
    matrix_ls = []
    for seq_sple, comp_path in tqdm(zip(sample_info['seq_sple'],
            sample_info['dna_sample_comp_path'])):
        sample_ls = []
        parser = SampleCompParser(comp_path, ncbi_tax=ncbi)
        for taxid in panel_orgs['taxid']:
            nr = parser.get_taxid_nr(taxid, normalizer=1)
            sample_ls.append(nr)
        matrix_ls.append(np.array(sample_ls))
    return np.array(matrix_ls)


ncbi = NcbiTaxonomy()

panel_orgs = pd.read_csv('org_taxids_uti_all_with_new.txt', sep='\t')
samples_arup = pd.read_csv('arup_sample_info.txt', sep='\t')
samples_syn = pd.read_csv('synergy_clinical_sample_info.txt', sep='\t')

cov_arup = pd.read_csv('arup_coverages.csv', index_col=2)
cov_arup_fp = pd.read_csv('arup_coverages_fungpar.csv', index_col=2)
cov_syn = pd.read_csv('synergy_coverages.csv', index_col=2)
cov_syn_fp = pd.read_csv('synergy_coverages_fungpar.csv', index_col=2)

merged_arup = pd.concat([cov_arup, cov_arup_fp])
merged_syn = pd.concat([cov_syn, cov_syn_fp])

quant_arup = get_quant(samples_arup, panel_orgs, merged_arup)
quant_arup_df = pd.DataFrame(data=quant_arup, columns=panel_orgs['name'],
    index=samples_arup['seq_sple'])
quant_arup_df.to_csv('sample_comp_abundance_on_panel_arup.csv')

quant_syn = get_quant(samples_syn, panel_orgs, merged_syn)
quant_syn_df = pd.DataFrame(data=quant_syn, columns=panel_orgs['name'],
    index=samples_syn['seq_sple'])
quant_syn_df.to_csv('sample_comp_abundance_on_panel_syn.csv')
