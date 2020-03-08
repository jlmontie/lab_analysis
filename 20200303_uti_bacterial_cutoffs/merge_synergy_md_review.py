import pandas as pd
import numpy as np

coverages = pd.read_csv('synergy_coverages.csv',
    dtype={'accession': str, 'seq_id': str, 'batch_id': str, 'seq_sple': str,
        'organism': str, 'taxid': str},
    index_col=0)
dump_orgs = pd.read_excel('explify-synergy-export-200303.xlsx',
    sheet_name='organisms', dtype={'Accession': str})
dump_orgs = dump_orgs.rename(columns={'Accession': 'accession',
    'Organism Name': 'organism', 'Review': 'review', 'Batch ID': 'batch_id',
    'Reporting ID': 'reporting_id'})
merged = coverages.merge(
    dump_orgs[['batch_id', 'accession', 'organism', 'reporting_id', 'review']],
    on=['batch_id', 'accession', 'organism'], how='left'
)
print(merged['review'].value_counts())
merged['review'] = merged['review'].replace({
    'DETECTED': 'positive', 'REJECTED': 'rejected',
    'INCONCLUSIVE': 'inconclusive'
})
print(merged['review'].value_counts())
merged['above_cutoff'] = True
merged.loc[merged['reporting_id'].isna(), 'above_cutoff'] = False
merged_filtered = merged.drop_duplicates(subset=['accession', 'batch_id', 'organism'])
merged_filtered = merged_filtered.sort_values(['batch_id'])
print(f"merged.shape before: {merged_filtered.shape}")
merged_filtered = merged_filtered.drop_duplicates(subset=['accession', 'organism'], keep='last')
print(f"merged.shape after: {merged_filtered.shape}")
print(f"merged no na reporting_id: {np.sum(~merged_filtered['reporting_id'].isna())}")
print(f"merged no na above_cutoff: {np.sum(merged_filtered['above_cutoff'])}")
print(f"merged columns: {merged_filtered.columns}")
merged_filtered.to_csv('synergy_coverages_with_review.csv', index=False)
