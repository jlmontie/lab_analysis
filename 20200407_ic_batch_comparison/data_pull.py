from collections import defaultdict

import pandas as pd

from data_extractor import DataExtractor


def get_accession_dict(df):
    accession_dict = defaultdict(list)
    for _, row in df.iterrows():
        accession_dict[row['batch_id']].append(row['accession'])
    return accession_dict


def get_data(project_dirs, batch_dict):
    df_ls = []
    batches_found = []
    batches_not_found = []
    for project_dir in project_dirs:
        extractor = DataExtractor(project_dir=project_dir,
            accession_dict=batch_dict, lib_type='rna')
        df_project_dir = extractor.collect_data()
        df_ls.append(df_project_dir)
        batches_found.extend(list(extractor.batches_found))
        batches_not_found.extend(list(extractor.batches_not_found))
    df = pd.concat(df_ls)
    # batches not found in one project directory may be found in another
    batches_not_found_all = list(
        set(batches_not_found).difference(batches_found))
    return df, batches_found, batches_not_found_all

project_dir_1 = '/srv/idbydna-group3/results/idbd_rnd_v2/'
project_dir_2 = '/srv/idbydna-group2/results/arup-resp-prod/'
project_dirs = [project_dir_1, project_dir_2]

batch_info = pd.read_csv('ic_batch_info.csv')
batch_new = batch_info[batch_info['ic_batch'] == 'new']
batch_new_dict = get_accession_dict(batch_new)
batch_old = batch_info[batch_info['ic_batch'] == 'old']
batch_old_dict = get_accession_dict(batch_old)

######## Old batches ########
df_old, batches_found_old, batches_not_found_old = get_data(
    project_dirs, batch_old_dict
)
df_old.to_csv('old_ic_batch.csv', index=False)
print(df_old)
print(f"Batches found for old:\n{batches_found_old}\n")
print(f"Batches not found for old:\n{batches_not_found_old}")

######## New batches ########
df_new, batches_found_new, batches_not_found_new = get_data(
    project_dirs, batch_new_dict
)
df_new.to_csv('new_ic_batch.csv', index=False)
print(df_new)
print(f"Batches found for new:\n{batches_found_new}\n")
print(f"Batches not found for new:\n{batches_not_found_new}")
