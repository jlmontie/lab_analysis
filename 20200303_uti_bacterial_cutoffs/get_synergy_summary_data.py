import json
import gzip
from collections import defaultdict

import pandas as pd
from tqdm import tqdm

info = pd.read_csv('synergy_clinical_sample_info.txt', sep='\t')
rna_data = defaultdict(list)
dna_data = defaultdict(list)
for idx, row in tqdm(info.iterrows()):
    rna_path = row['rna_bac_summary_path']
    with gzip.open(rna_path) as infile:
        for line in infile:
            obj = json.loads(line)
            rna_data['seq_id'].append(row['seq_id'])
            rna_data['batch_id'].append(row['batch_id'])
            rna_data['seq_sple'].append(row['seq_sple'])
            rna_data['accession'].append(row['accession'])
            rna_data['taxid'].append(obj['taxid'])
            rna_data['organism'].append(obj['name'])
            rna_data['rna_coverage'].append(obj['coverage'])
            rna_data['absolute_quant'].append(obj['absolute_quant'])
    dna_path = row['dna_bac_summary_path']
    with gzip.open(dna_path) as infile2:
        for line in infile2:
            obj = json.loads(line)
            dna_data['seq_id'].append(row['seq_id'])
            dna_data['batch_id'].append(row['batch_id'])
            dna_data['seq_sple'].append(row['seq_sple'])
            dna_data['accession'].append(row['accession'])
            dna_data['taxid'].append(obj['taxid'])
            dna_data['organism'].append(obj['name'])
            dna_data['dna_coverage'].append(obj['coverage'])
            dna_data['core_coverage'].append(obj['core_coverage'])
rna_df = pd.DataFrame(rna_data)
dna_df = pd.DataFrame(dna_data)
merged = rna_df.merge(dna_df,
    on=['seq_id', 'batch_id', 'seq_sple', 'accession', 'taxid', 'organism'],
    how='outer')
print(f"rna_df.shape: {rna_df.shape}")
print(f"dna_df.shape: {dna_df.shape}")
print(f"merged.shape: {merged.shape}")
merged.to_csv('synergy_coverages.csv', index=False)
