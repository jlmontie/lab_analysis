from collections import defaultdict

import pandas as pd

from data_extractor import DataExtractor


def get_accession_dict(df):
    accession_dict = defaultdict(list)
    for _, row in df.iterrows():
        accession_dict[row['batch_id']].append(row['accession'])
    return accession_dict


project_dir_1 = '/srv/idbydna-group3/results/idbd_rnd_v2/'
project_dir_2 = '/srv/idbydna-group2/results/arup-resp-prod/'
project_dirs = [project_dir_1, project_dir_2]

batch_info = pd.read_csv('ic_batch_info.csv')
batch_new = batch_info[batch_info['ic_batch'] == 'new']
batch_new_dict = get_accession_dict(batch_new)
batch_old = batch_info[batch_info['ic_batch'] == 'old']
batch_old_dict = get_accession_dict(batch_old)

df_old_ls = []
for project_dir in project_dirs:
    extractor_old = DataExtractor(project_dir=project_dir,
        accession_dict=batch_old_dict, lib_type='rna')
    df_old_ls.append(extractor_old.collect_data())
df_old = pd.concat(df_old_ls)
print(df_old)
