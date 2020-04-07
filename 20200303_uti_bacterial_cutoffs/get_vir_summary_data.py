import json
import gzip
import argparse
from collections import defaultdict

import pandas as pd
from tqdm import tqdm


def main(input, output):
    info = pd.read_csv(input, sep='\t')
    rna_data = defaultdict(list)
    dna_data = defaultdict(list)
    for idx, row in tqdm(info.iterrows()):
        rna_path = row['rna_vir_summary_path']
        if rna_path.endswith('.gz'):
            infile = gzip.open(rna_path)
        else:
            infile = open(rna_path)
        for line in infile:
            obj = json.loads(line)
            rna_data['seq_id'].append(row['seq_id'])
            rna_data['batch_id'].append(row['batch_id'])
            rna_data['seq_sple'].append(row['seq_sple'])
            rna_data['accession'].append(row['accession'])
            rna_data['taxid'].append(obj['taxid'])
            rna_data['organism'].append(obj['name'])
            rna_data['rna_coverage'].append(obj['coverage'])
        dna_path = row['dna_vir_summary_path']
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
    rna_df = pd.DataFrame(rna_data)
    dna_df = pd.DataFrame(dna_data)
    merged = rna_df.merge(dna_df,
        on=['seq_id', 'batch_id', 'seq_sple', 'accession', 'taxid', 'organism'],
        how='outer')
    print(f"rna_df.shape: {rna_df.shape}")
    print(f"dna_df.shape: {dna_df.shape}")
    print(f"merged.shape: {merged.shape}")
    merged.to_csv(output, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    parser.add_argument('output')
    args = parser.parse_args()
    main(args.input, args.output)