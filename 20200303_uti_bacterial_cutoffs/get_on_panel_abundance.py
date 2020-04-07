import sys
from collections import defaultdict

import pandas as pd
import numpy as np
from tqdm import tqdm

from idbd_bio_utils import NcbiTaxonomy
sys.path.insert(0, '../scripts')
from sample_composition_utils import SampleCompParser


def make_sample_dict(coverages):
    coverages_dict = defaultdict(dict)
    for idx, row in coverages.iterrows():
        coverages_dict[row['seq_sple']].update(
            {row['taxid']: row['absolute_quant']}
        )
    return coverages_dict


def get_quant(sample_info, panel_orgs, coverages):
    coverages_dict = make_sample_dict(coverages)
    matrix_abs = []
    matrix_comp = []
    for seq_sple, comp_path in tqdm(zip(sample_info['seq_sple'],
            sample_info['dna_sample_comp_path'])):
        sample_abs = []
        sample_comp = []
        parser = SampleCompParser(comp_path, ncbi_tax=ncbi)
        for taxid in panel_orgs['taxid']:
            if taxid in coverages_dict[seq_sple]:
                quant_abs = coverages_dict[seq_sple][taxid]
                if np.isnan(quant_abs):
                    quant_abs = 0
                nr = parser.get_taxid_nr(taxid, normalizer=1)
            else:
                quant_abs = 0
                nr = 0
            sample_abs.append(quant_abs)
            sample_comp.append(nr)
        matrix_abs.append(np.array(sample_abs))
        matrix_comp.append(np.array(sample_comp))
    return np.array(matrix_abs), np.array(matrix_comp)


ncbi = NcbiTaxonomy()

panel_orgs = pd.read_csv('org_taxids_uti_all_with_new.txt', sep='\t')
samples_arup = pd.read_csv('arup_sample_info.txt', sep='\t')
samples_syn = pd.read_csv('synergy_clinical_sample_info.txt', sep='\t')

cov_arup = pd.read_csv('arup_coverages_with_review.csv', index_col=2)
cov_arup_fp = pd.read_csv('arup_coverages_fungpar_with_review.csv', index_col=2)
cov_syn = pd.read_csv('synergy_coverages_with_review.csv', index_col=2)
cov_syn_fp = pd.read_csv('synergy_coverages_fungpar_with_review.csv', index_col=2)

merged_arup = pd.concat([cov_arup, cov_arup_fp])
merged_arup = merged_arup[merged_arup['above_cutoff']]
merged_syn = pd.concat([cov_syn, cov_syn_fp])
merged_syn = merged_syn[merged_syn['above_cutoff']]

quant_arup_abs, quant_arup_comp = get_quant(samples_arup, panel_orgs,
    merged_arup)
quant_arup_df_abs = pd.DataFrame(data=quant_arup_abs, columns=panel_orgs['name'],
    index=samples_arup['seq_sple'])
quant_arup_df_abs.to_csv('absolute_quant_on_panel_arup.csv')
quant_arup_df_comp = pd.DataFrame(data=quant_arup_comp, columns=panel_orgs['name'],
    index=samples_arup['seq_sple'])
quant_arup_df_comp.to_csv('sample_comp_abundance_on_panel_arup.csv')

quant_syn_abs, quant_syn_comp = get_quant(samples_syn, panel_orgs,
    merged_syn)
quant_syn_df_abs = pd.DataFrame(data=quant_syn_abs, columns=panel_orgs['name'],
    index=samples_syn['seq_sple'])
quant_syn_df_abs.to_csv('absolute_quant_on_panel_syn.csv')
quant_syn_df_comp = pd.DataFrame(data=quant_syn_comp, columns=panel_orgs['name'],
    index=samples_syn['seq_sple'])
quant_syn_df_comp.to_csv('sample_comp_abundance_on_panel_syn.csv')
