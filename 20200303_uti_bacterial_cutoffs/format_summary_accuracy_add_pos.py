import pandas as pd
import numpy as np

report = pd.read_csv('Summary_Accuracy_200226.csv')
print(report.columns)
add_pos = report['Additional Positives'].str.split(', ', expand=True)
joined = report.join(add_pos)
print(joined)
melt = joined.melt(id_vars=[
    'Sample ID', 'Accession number', 'Batch', 'Match?'
    ],
    value_vars=['ARUP results', 0, 1, 2, 3, 4, 5],
    value_name='organism'
)
melt.loc[melt['variable'].isin(
    [0, 1, 2, 3, 4, 5]), 'Match?'] = 'Additional Pos.'
melt = melt.drop(columns=['variable'])
melt = melt[~melt['organism'].isin([None, np.nan])]
melt = melt.rename(columns={'Accession number': 'accession',
                            'Match?': 'review',
                            'Batch': 'batch_id',
                            'Sample ID': 'sample_id'})
melt = melt[melt['organism'].str.contains(' ')]
melt = melt[~(melt['organism'] == 'No growth')]
print(melt)
melt.to_csv('Summary_Accuracy_200226_formatted.csv', index=False)
