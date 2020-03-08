import pandas as pd
from tqdm import tqdm

coverages = pd.read_csv('coverages_all.csv',
    dtype={'accession': str, 'seq_sple': str, 'batch_id': str, 'seq_sple': str,
        'organism': str, 'taxid': int, 'review': str, 'above_cutoff': str,
        'source': str})
coverages = coverages.sort_values('organism')
with open('org_taxids_all.txt', 'w') as outfile:
    outfile.write(f"name\ttaxid\n")
    for taxid in tqdm(coverages['taxid'].unique()):
        df_tmp = coverages[coverages['taxid'] == taxid]
        try:
            most_common_name = df_tmp['organism'].mode().values[0]
        except:
            print(df_tmp['organism'])
            print(df_tmp['organism'].mode())
        outfile.write(f"{most_common_name}\t{taxid}\n")