import os
import glob
import json
import gzip

import pandas as pd
import numpy as np

from bin_sample_composition import bin_reads
from ncbi_taxonomy_utils import ncbi_taxonomy

ncbi = ncbi_taxonomy()


def parse_composition_file(composition_path, taxid):
    with open(composition_path) as comp_file:
        total_read_cnt = 0
        taxid_cnt = []
        for line in comp_file:
            data = line.strip().split('\t')
            data = [int(item) for item in data]
            total_read_cnt += data[1]
            if isinstance(taxid, list):
                if data[0] in taxid:
                    taxid_cnt.append((data[0], data[1]))
            else:
                if data[0] == taxid:
                    taxid_cnt.append((data[0], data[1]))
    return_dict = {
        'Total Read Count': total_read_cnt,
        'Taxid Read Count': taxid_cnt
    }
    return return_dict


def get_lib_dict_ls(results_dir=None, fqo=None):
    fqo_hi_ic = pd.read_csv(fqo)
    fqo_hi_ic['Run Date (YYYY-MM-DD)'] = pd.to_datetime(fqo_hi_ic['Run Date (YYYY-MM-DD)'])
    results_dir = results_dir
    batch_paths = glob.glob(os.path.join(results_dir, 'batch', '*'))
    lib_dict_ls = []
    for batch_path in batch_paths:
        with open(batch_path) as batch_file:
            batch = json.load(batch_file)
            batch_id = batch['batch']['libBatchId']
        for lib in batch['libraries']:
            accession = lib['bioSple']
            seq_sple = lib['seqSple']
            lib_type = lib['libType']
            if not any(fqo_hi_ic['Seq Sple'] == seq_sple.lower()):
                continue
            run_date = fqo_hi_ic.loc[fqo_hi_ic['Seq Sple'] == seq_sple.lower(), 'Run Date (YYYY-MM-DD)'].values[0]
            sample_name = lib['spleName']
            total_reads = lib['qualityFilterInfo']['readsOut']
            tax_paths = lib['diagnosticOutput']
            vir_paths = [path for path in tax_paths if 'dna.viral.dxsm.out.summary.gz' in path]
            t7_read_cnt = 0
            pr772_read_cnt = 0
            if vir_paths:
                vir_path = vir_paths[0]
                with gzip.open(os.path.join(results_dir, vir_path)) as vir_file:
                    for line in vir_file:
                        obj = json.loads(line.strip())
                        if obj['reporting_id'] == '26706_10760':
                            t7_read_cnt = obj['read_count']
                        elif obj['reporting_id'] == '26648_261665':
                            pr772_read_cnt = obj['read_count']
            t7_norm = 10e6 * t7_read_cnt / total_reads
            pr772_norm = 10e6 * pr772_read_cnt / total_reads
            # composition file data
            composition_path_ls = glob.glob(os.path.join(results_dir, 'tax', seq_sple + '*dna.sample_composition.out'))
            if not composition_path_ls:
                composition_data = {
                    'Total Read Count': 1e-6,
                    'Taxid Read Count': [(10760, 0), (261665, 0)]
                }
                org_composition = {
                    "Human": 0,
                    "Bacteria": 0,
                    "Virus": 0,
                    "Parasite": 0,
                    "Fungus": 0,
                    "Unclassified": 0
                }
                composition_file = None
            else:
                composition_path = composition_path_ls[0]
                composition_file = os.path.basename(composition_path)
                composition_data = parse_composition_file(composition_path, [10760, 261665])
                org_composition = bin_reads(composition_path, ncbi_class=ncbi, quantification='relative', ctrl_taxids=[10760, 261665])
                # for key in org_composition:
                #     org_composition[key] = 10e6 * org_composition[key] / composition_data['Total Read Count']
            if not composition_data['Taxid Read Count']:
                composition_data['Taxid Read Count'] = [(10760, 0), (261665, 0)]
            elif len(composition_data['Taxid Read Count']) == 1:
                if composition_data['Taxid Read Count'][0][0] == 10760:
                    composition_data['Taxid Read Count'].append((261665, 0))
                elif composition_data['Taxid Read Count'][0][0] == 261665:
                    composition_data['Taxid Read Count'].append((10760, 0))
            t7_raw_comp = [item[1] for item in composition_data['Taxid Read Count'] if item[0] == 10760][0]
            pr772_raw_comp = [item[1] for item in composition_data['Taxid Read Count'] if item[0] == 261665][0]
            t7_norm_comp = 10e6 * t7_raw_comp / composition_data['Total Read Count']
            pr772_norm_comp = 10e6 * pr772_raw_comp / composition_data['Total Read Count']
            lib_dict = {
                'Accession': accession,
                'Seq Sple': seq_sple,
                'Composition File': composition_file,
                'Batch ID': batch_id,
                'Run Date': run_date,
                'Sample Name': sample_name,
                'Library Type': lib_type,
                'Total Reads': total_reads,
                'T7 Raw Reads': t7_read_cnt,
                'PR772 Raw Reads': pr772_read_cnt,
                'T7 Normalized Reads': t7_norm,
                'PR772 Normalized Reads': pr772_norm,
                'T7 + PR772 NR': t7_norm + pr772_norm,
                'Log10 T7 + PR772 NR': np.log10(t7_norm + pr772_norm),
                'Total Reads (Composition File)': composition_data['Total Read Count'],
                'T7 Raw Reads (Composition File)': t7_raw_comp,
                'T4 Raw Reads (Composition File)': pr772_raw_comp,
                'T7 Normalized (Composition File)': t7_norm_comp,
                'T4 Normalized (Composition File)': pr772_norm_comp,
                'T7 + PR772 NR (Composition File)': t7_norm_comp + pr772_norm_comp,
                'Log10 T7 + PR772 NR (Composition File)': np.log10(t7_norm_comp + pr772_norm_comp),
                'Summary - Composition (log10)': np.log10(t7_norm + pr772_norm) - np.log10(t7_norm_comp + pr772_norm_comp),
                "Human": org_composition['Human'],
                "Bacteria": org_composition['Bacteria'],
                "Virus": org_composition['Virus'],
                "Parasite": org_composition['Parasite'],
                "Fungus": org_composition['Fungus'],
                "Unclassified": org_composition['Unclassified']
            }
            lib_dict_ls.append(lib_dict)
    return lib_dict_ls

# Diluted Hi IC
hi_ic_fqo = 'FastQataloguer_HighIC_2020-01-22.csv'
hi_ic_results_dir = '/data/analysis_group1/idbd_rnd/results/200120_NB551543_0200_AH2N2NAFX2'
hi_ic_lib_dict_ls = get_lib_dict_ls(results_dir=hi_ic_results_dir, fqo=hi_ic_fqo)
hi_ic_df = pd.DataFrame(hi_ic_lib_dict_ls)
hi_ic_df.to_csv('high_ic_results.csv', index=False)

# Low IC with a new negative control
# First batch
new_low_ic_fqo = 'FastQataloguer_191210b02_2020-01-25.csv'
new_low_ic_results_dir = '/data/analysis_group1/idbd_rnd/results/191211_NB551543_0175_AHVYMFAFXY'
new_low_ic_lib_dict_ls = get_lib_dict_ls(results_dir=new_low_ic_results_dir, fqo=new_low_ic_fqo)
new_low_ic_df = pd.DataFrame(new_low_ic_lib_dict_ls)
new_low_ic_df.to_csv('new_low_ic_results.csv', index=False)
# Second batch
new_low_ic_fqo_2 = 'FastQataloguer_191211b02_2020-01-26.csv'
new_low_ic_results_dir_2 = '/data/analysis_group1/idbd_rnd/results/191213_NB551702_0091_AHW3MKAFXY'
new_low_ic_lib_dict_ls_2 = get_lib_dict_ls(results_dir=new_low_ic_results_dir_2, fqo=new_low_ic_fqo_2)
new_low_ic_df_2 = pd.DataFrame(new_low_ic_lib_dict_ls_2)
new_low_ic_df_2.to_csv('new_low_ic_results_2.csv', index=False)

# Subsequent Low IC
fqo_lo_ic = pd.read_csv('FastQataloguer_LowIC_2020-01-22.csv')
fqo_lo_ic['Run Date (YYYY-MM-DD)'] = pd.to_datetime(fqo_lo_ic['Run Date (YYYY-MM-DD)'])
low_ic_comp_dict_ls = []
for low_ic_comp in os.listdir('low_ic_compositions'):
    seqid = '-'.join(low_ic_comp.split('-')[:3])
    seq_sple = low_ic_comp.split('.')[0]
    run_date = fqo_lo_ic.loc[(fqo_lo_ic['Sequencing Id'] == seqid) & (fqo_lo_ic['Accession'] == 'controlneg'),
                             'Run Date (YYYY-MM-DD)'].values[0]
    accession = fqo_lo_ic.loc[(fqo_lo_ic['Sequencing Id'] == seqid) & (fqo_lo_ic['Accession'] == 'controlneg'),
                              'Accession'].values[0]
    batch_id = fqo_lo_ic.loc[(fqo_lo_ic['Sequencing Id'] == seqid) & (fqo_lo_ic['Accession'] == 'controlneg'),
                             'Lib Batch Id'].str.upper().values[0]
    low_ic_comp_data = parse_composition_file(os.path.join('low_ic_compositions', low_ic_comp), [10760, 261665])
    low_ic_t7_raw_comp = [item[1] for item in low_ic_comp_data['Taxid Read Count'] if item[0] == 10760][0]
    low_ic_pr772_raw_comp = [item[1] for item in low_ic_comp_data['Taxid Read Count'] if item[0] == 261665][0]
    low_ic_t7_norm_comp = 10e6 * low_ic_t7_raw_comp / low_ic_comp_data['Total Read Count']
    low_ic_pr772_norm_comp = 10e6 * low_ic_pr772_raw_comp / low_ic_comp_data['Total Read Count']
    org_composition = bin_reads(os.path.join('low_ic_compositions', low_ic_comp),
                                ncbi_class=ncbi, quantification='relative',
                                ctrl_taxids=[10760, 261665])
    # for key in org_composition:
    #     org_composition[key] = 10e6 * org_composition[key] / low_ic_comp_data['Total Read Count']
    low_ic_comp_dict = {
        'Accession': accession,
        'Seq Sple': seq_sple,
        'Composition File': low_ic_comp,
        'Batch ID': batch_id,
        'Run Date': run_date,
        'T7 Raw Reads (Composition File)': low_ic_t7_raw_comp,
        'T4 Raw Reads (Composition File)': low_ic_pr772_raw_comp,
        'T7 Normalized (Composition File)': low_ic_t7_norm_comp,
        'T4 Normalized (Composition File)': low_ic_pr772_norm_comp,
        'T7 + PR772 NR (Composition File)': low_ic_t7_norm_comp + low_ic_pr772_norm_comp,
        'Log10 T7 + PR772 NR (Composition File)': np.log10(low_ic_t7_norm_comp + low_ic_pr772_norm_comp),
        "Human": org_composition['Human'],
        "Bacteria": org_composition['Bacteria'],
        "Virus": org_composition['Virus'],
        "Parasite": org_composition['Parasite'],
        "Fungus": org_composition['Fungus'],
        "Unclassified": org_composition['Unclassified']
    }
    low_ic_comp_dict_ls.append(low_ic_comp_dict)

low_ic_df = pd.DataFrame(low_ic_comp_dict_ls)
low_ic_df.to_csv('low_ic_results.csv', index=False)
