from collections import defaultdict

import pandas as pd

from data_extractor import DataExtractor


def get_accession_dict(df):
    accession_dict = defaultdict(list)
    for _, row in df.iterrows():
        accession_dict[row['batch_id']].append(row['accession'])
    return accession_dict


project_dir = '/srv/idbydna-group3/results/idbd_rnd_v2/'

batch_info = pd.read_csv('ic_batch_info.csv')
batch_new = batch_info[batch_info['ic_batch'] == 'new']
batch_new_dict = get_accession_dict(batch_new)
batch_old = batch_info[batch_info['ic_batch'] == 'old']
batch_old_dict = get_accession_dict(batch_old)

extractor_old = DataExtractor(project_dir=project_dir,
    accession_dict=batch_old_dict, lib_type='rna')
df_old = extractor_old.collect_data()
print(df_old)
