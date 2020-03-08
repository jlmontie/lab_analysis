import pandas as pd
from tqdm import tqdm

coverages = pd.read_csv('coverages_all.csv',
    dtype={'accession': str, 'seq_sple': str, 'batch_id': str, 'seq_sple': str,
        'organism': str, 'taxid': int, 'review': str, 'above_cutoff': str,
        'source': str})
for taxid in tqdm(coverages['taxid'].unique()):
    save_df = coverages[coverages['taxid'] == taxid]
    save_df.to_csv(f"data/{taxid}.csv", index=False)