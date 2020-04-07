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
        rna_path = row['rna_fungpar_summary_path']
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
    rna_df = pd.DataFrame(rna_data)
    print(f"rna_df.shape: {rna_df.shape}")
    rna_df.sort_values('rna_coverage', ascending=False).to_csv(output,
                                                               index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    parser.add_argument('output')
    args = parser.parse_args()
    main(args.input, args.output)