import os
import glob
import json

import pandas as pd

run_paths = []
with open('run_file_arup.txt') as infile:
    for line in infile:
        run_paths.append(line.strip())

run_path_ls = []
run_id_ls = []
batch_id_ls = []
accession_ls = []
seq_sple_ls = []
study_ls = []
spletype_ls = []
bac_dna_summary_path_ls = []
bac_rna_summary_path_ls = []
bac_amr_summary_path_ls = []
run_date_ls = []
for run_path in run_paths:
    batch_dir = os.path.join(run_path, 'batch')
    batch_paths = glob.glob(os.path.join(batch_dir, '*'))
    for batch_path in batch_paths:
        with open(batch_path) as batch_infile:
            batch = json.load(batch_infile)
        for lib in batch['libraries']:
            if lib['bioSple'] in ['controlBlk', 'controlPos', 'controlNeg']:
                continue
            run_path_ls.append(run_path)
            run_id_ls.append(batch['sequencing']['sequencingId'])
            batch_id_ls.append(batch['batch']['libBatchId'])
            accession_ls.append(lib['bioSple'])
            seq_sple_ls.append(lib['seqSple'])
            study_ls.append(lib['study'])
            spletype_ls.append(lib['spleType'])
            bac_rna = [os.path.join(run_path, item) for item in lib['rnaDiagnosticOutput'] if 'rna.bacterial.dxsm' in item]
            if bac_rna:
                bac_rna_summary_path_ls.extend(bac_rna)
            else:
                bac_rna_summary_path_ls.append(None)
            bac_dna = [os.path.join(run_path, item) for item in lib['dnaDiagnosticOutput'] if 'dna.bacterial.dxsm' in item]
            if bac_dna:
                bac_dna_summary_path_ls.extend(bac_dna)
            else:
                bac_dna_summary_path_ls.append(None)
            bac_amr = [os.path.join(run_path, item) for item in lib['dnaDiagnosticOutput'] if 'dna.bacterial_amr.dxsm' in item]
            if bac_amr:
                bac_amr_summary_path_ls.extend(bac_amr)
            else:
                bac_amr_summary_path_ls.append(None)
            run_date_ls.append(batch['analysis']['timeCompleted'][:10])

assert all((len(run_path_ls) == len(run_id_ls), len(run_path_ls) == len(accession_ls), len(run_path_ls) == len(seq_sple_ls),
    len(run_path_ls) == len(bac_rna_summary_path_ls), len(run_path_ls) == len(bac_dna_summary_path_ls),
    len(run_path_ls) == len(run_date_ls), len(run_path_ls) == len(spletype_ls), len(run_path_ls) == len(study_ls),
    len(run_path_ls) == len(bac_amr_summary_path_ls)))

df = pd.DataFrame({
    'run_dir': run_path_ls,
    'seq_id': run_id_ls,
    'batch_id': batch_id_ls,
    'seq_sple': seq_sple_ls,
    'accession': accession_ls,
    'study': study_ls,
    'sample_type': spletype_ls,
    'rna_bac_summary_path': bac_rna_summary_path_ls,
    'dna_bac_summary_path': bac_dna_summary_path_ls,
    'dna_amr_summary_path': bac_amr_summary_path_ls,
    'run_date': run_date_ls
})
df = df.sort_values('run_date')
df.to_csv('arup_sample_info.txt', index=False, sep='\t')
